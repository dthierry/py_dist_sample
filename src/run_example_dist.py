# !/usr/bin/env python
# -*- coding:utf:8 -*-

from model.dmod import mod
from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import TerminationCondition, SolverFactory, SolverStatus
import platform
import sys

__author__ = "David Thierry"
"""
Set up the model.
Find initial condition.
Solve dynamic model.
Note: we need ipopt! (tested with hsl ma57)
"""


def main():
    if platform.python_version().startswith('2'):
        sys.stdout.write("Python 2.x detected! Python 3.x+ required!")
        sys.exit()
    #: Discretize.
    t = TransformationFactory("dae.collocation")
    t.apply_to(mod, nfe=1, ncp=1)

    #: Set up bounds.
    vml = 0.5 * ((1 / 2288) * 0.2685 ** (1 + (1 - 100 / 512.4) ** 0.2453)) + \
          0.5 * ((1 / 1235) * 0.27136 ** (1 + (1 - 100 / 536.4) ** 0.24))
    vmu = 0.5 * ((1 / 2288) * 0.2685 ** (1 + (1 - (512.4 - 1E-08) / 512.4) ** 0.2453)) + \
          0.5 * ((1 / 1235) * 0.27136 ** (1 + (1 - 512.4 / 536.4) ** 0.24))

    mod.M[:, :].setlb(-1E-08)
    mod.L[:, :].setlb(-1E-08)
    mod.V[:, :].setlb(-1E-08)

    mod.Mv[:, :].setlb(0.155 - 1E-08)
    mod.Mv1[:].setlb(8.5 - 1E-08)
    mod.Mvn[:].setlb(0.17 - 1E-08)
    #
    mod.T[:, :].setlb(100.0)
    mod.T[:, :].setub(512.4 - 1E-08)

    mod.y[:, :].setlb(-1E-11)
    mod.x[:, :].setlb(-1E-11)

    mod.y[:, :].setub(1.0 + 1E-11)
    mod.x[:, :].setub(1.0 + 1E-11)
    mod.Vm[:, :].setlb(vml)
    mod.Vm[:, :].setub(vmu)
    ####
    #: Solve steady state
    ####
    #: Deactivate discretization equations.
    mod.xdot_disc_eq.deactivate()
    mod.Mdot_disc_eq.deactivate()
    mod.xdot[:, :].fix(0.0)
    mod.Mdot[:, :].fix(0.0)
    #: Deactivate initial conditions.
    mod.M_icc.deactivate()
    mod.x_icc.deactivate()
    #: Try with the rest of the constraints
    for i in mod.component_objects(Constraint):
        s = i.index_set()
        if s.dimen > 2:
            s.set_tuple()
            for ss in s.set_tuple():
                if ss is mod.t:
                    ss[0, :].deactivate()
        else:
            if s is mod.t:
                try:
                    i[0].deactivate()
                except KeyError:
                    pass

    #: Set up a dummy objective.
    mod.o = Objective(expr=(mod.u2[mod.t.last()] - value(mod.Qr[mod.t.last()])) ** 2)
    #: Set up ipopt.
    ip = SolverFactory("/home/dav0/apps/ipopt-build-vanilla/bin/ipopt")  #: Ipopt's executable binary goes here.
    ip.options["halt_on_ampl_error"] = "yes"
    ip.options["linear_solver"] = "ma57"
    res = ip.solve(mod, tee=True, symbolic_solver_labels=True)
    if res.solver.status == SolverStatus.ok and res.solver.termination_condition == TerminationCondition.optimal:
        print("Good")
    else:
        print("error!")
        sys.exit()

    #: Print to files.
    mod.pprint(filename="model_steady.pprint")
    mod.display(filename="model_steady.display")
    #: Done solving.

    #: Set up dynamic problem initial values
    uval = value(mod.Qr[1])
    for v in mod.Qr.itervalues():
        v.set_value(uval)
        # v.fix()
    for v in mod.u2.itervalues():
        v.set_value(uval)
        v.fix()

    #: Set initial conditions
    for tray in mod.tray:
        mod.x_ic[tray].set_value(value(mod.x[1, tray]))
        mod.M_ic[tray].set_value(value(mod.M[1, tray]))

    mod.xdot_disc_eq.activate()
    mod.Mdot_disc_eq.activate()
    mod.xdot[:, :].unfix()
    mod.Mdot[:, :].unfix()
    #: Activate initial conditions
    mod.M_icc.activate()
    mod.x_icc.activate()
    #: Try with the rest of the constraints
    for i in mod.component_objects(Constraint):
        s = i.index_set()
        if s.dimen > 2:
            s.set_tuple()
            for ss in s.set_tuple():
                if ss is mod.t:
                    ss[0, :].activate()
        else:
            if s is mod.t:
                try:
                    i[0].activate()
                except KeyError:
                    pass
    mod.o.deactivate()  #: Forget about the dummy objective.
    ip.solve(mod, tee=True, symbolic_solver_labels=False)
    mod.display(filename="dynamic.display")
    ip.options["tol"] = 1E-06
    res = ip.solve(mod, tee=True, symbolic_solver_labels=False)

    if res.solver.status == SolverStatus.ok and res.solver.termination_condition == TerminationCondition.optimal:
        print("Good. Check the dynamic.display file!")
    else:
        print("error!")
        sys.exit()


if __name__ == '__main__':
    main()