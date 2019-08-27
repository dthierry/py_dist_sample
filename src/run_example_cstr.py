#!/usr/bin/env python
# -*- coding: utf-8 -*-

from cstr.cstr_ import mod
from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import TerminationCondition, SolverFactory, SolverStatus
import platform
import sys

__author__ = "David Thierry"


def main():
    if platform.python_version().startswith('2'):
        sys.stdout.write("Python 2.x detected! Python 3.x+ required!")
        sys.exit()

    t = TransformationFactory("dae.collocation")
    t.apply_to(mod, nfe=10, ncp=3)

    #: Set-up bounds
    #:

    #: Solve steady state
    ###
    mod.Cadot_disc_eq.deactivate()
    mod.Tdot_disc_eq.deactivate()
    mod.Tjdot_disc_eq.deactivate()

    mod.Cadot[:, :].fix(0.0)
    mod.Tdot[:, :].fix(0.0)
    mod.Tjdot[:, :].fix(0.0)
    #: Deactivate initial conditions
    mod.Ca_icc.deactivate()
    mod.T_icc.deactivate()
    mod.Tj_icc.deactivate()

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

    mod.o = Objective(expr=1, sense=minimize)
    ip = SolverFactory("ipopt")  #: Ipopt's executable binary goes here.
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
    uval = value(mod.Tjinb[1])
    for v in mod.Tjinb.itervalues():
        v.set_value(uval)
        # v.fix()

    #: Set initial conditions
    mod.Ca_ic[0].set_value(value(mod.Ca[1, 0]))
    mod.T_ic[0].set_value(value(mod.T[1, 0]))
    mod.Tj_ic[0].set_value(value(mod.Tj[1, 0]))

    mod.Cadot_disc_eq.activate()
    mod.Tdot_disc_eq.activate()
    mod.Tjdot_disc_eq.activate()
    mod.Cadot[:, :].unfix()
    mod.Tdot[:, :].unfix()
    mod.Tjdot[:, :].unfix()
    #: Activate initial conditions
    mod.Ca_icc.activate()
    mod.T_icc.activate()
    mod.Tj_icc.activate()
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


if __name__ == "__main__":
    main()
