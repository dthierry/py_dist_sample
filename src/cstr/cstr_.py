#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.core.base import ConcreteModel
from pyomo.dae import *
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.core.base.numvalue import value
from pyomo.core.base import Constraint, Set, Param, Var, Suffix
from pyomo.core.kernel import exp

__author__ = "David Thierry @dthierry"  #: March 2018


def _rule_k(m, i, n):
    if i == 0:
        return Constraint.Skip
    else:
        return m.k[i, n] == m.k0 * exp(-m.Er / m.T[i, n])


def _rule_ca(m, i, n):
    if i == 0:
        return Constraint.Skip
    else:
        rule = m.Cadot[i, n] == (m.F[i] / m.V) * (m.Cainb - m.Ca[i, n]) - 2 * m.k[i, n] * m.Ca[i, n] ** 2
        return rule


def _rule_t(m, i, n):
    if i == 0:
        return Constraint.Skip
    else:
        return m.Tdot[i, n] == (m.F[i] / m.V) * (m.Tinb - m.T[i, n]) + \
               2.0 * m.dH / (m.rho * m.Cp) * m.k[i, n] * m.Ca[i, n] ** 2 - \
               m.UA / (m.V * m.rho * m.Cp) * (m.T[i, n] - m.Tj[i, n])


def _rule_tj(m, i, n):
    if i == 0:
        return Constraint.Skip
    else:
        return m.Tjdot[i, n] == \
               (m.Fw[i] / m.Vw) * (m.Tjinb[i] - m.Tj[i, n]) + m.UA / (m.Vw * m.rhow * m.Cpw) * (m.T[i, n] - m.Tj[i, n])


def _rule_ca0(m, n):
    return m.Ca[0, n] == m.Ca_ic[n]


def _rule_t0(m, n):
    return m.T[0, n] == m.T_ic[n]


def _rule_tj0(m, n):
    return m.Tj[0, n] == m.Tj_ic[n]


#: type: (int, int, dict)
"""
    CSTR from Rodrigo's thesis
Returns:
    cstr_rodrigo_dae: The model itmod. Without discretization.
"""
#: if steady == True fallback to steady-state computation
mod = ConcreteModel()
#: mod.nfe_t = nfe_t  #:
#: mod.ncp_t = ncp_t
mod.discretized = False
ncstr = 1
mod.ncstr = Set(initialize=[i for i in range(0, ncstr)])

mod.t = ContinuousSet(bounds=(0, 1))

mod.Cainb = Param(default=1.0)
mod.Tinb = Param(default=275.0)
# mod.Tjinb = Param(default=250.0)

#: Our control var
mod.Tjinb = Var(mod.t, initialize=250)
mod.u1 = Param(mod.t, default=250, mutable=True)  #: We are making a sort-of port


def u1_rule(m, i):
    return m.Tjinb[i] == m.u1[i]


# mod.u1_cdummy = Constraint(mod.t, rule=lambda m, i: m.Tjinb[i] == mod.u1[i])
mod.u1_cdummy = Constraint(mod.t, rule=u1_rule)
#: u1 will contain the information from the NMPC problem. This is what drives the plant.
#: how about smth like nmpc_u1 or u1_nmpc

mod.V = Param(initialize=100)
mod.UA = Param(initialize=20000 * 60)
mod.rho = Param(initialize=1000)
mod.Cp = Param(initialize=4.2)
mod.Vw = Param(initialize=10)
mod.rhow = Param(initialize=1000)
mod.Cpw = Param(initialize=4.2)
mod.k0 = Param(initialize=4.11e13)
mod.E = Param(initialize=76534.704)
mod.R = Param(initialize=8.314472)
mod.Er = Param(initialize=lambda m: (value(mod.E) / value(mod.R)))
mod.dH = Param(initialize=596619.)

mod.F = Param(mod.t, mutable=True, default=1.2000000000000000E+02)
mod.Fw = Param(mod.t, mutable=True, default=3.0000000000000000E+01)

# States
mod.Ca = Var(mod.t, mod.ncstr, initialize=1.60659680385930765667001907104350E-02)
mod.T = Var(mod.t, mod.ncstr, initialize=3.92336059452774350120307644829154E+02)
mod.Tj = Var(mod.t, mod.ncstr, initialize=3.77995395658401662331016268581152E+02)

mod.k = Var(mod.t, mod.ncstr, initialize=4.70706140E+02)
mod.kdef = Constraint(mod.t, mod.ncstr)

#: These guys have to be zero at the steady-state (steady).
zero0 = dict.fromkeys(mod.t * mod.ncstr)
for key in zero0.keys():
    zero0[key] = 0.0
mod.steady = False
if mod.steady:
    mod.Cadot = zero0
    mod.Tdot = zero0
    mod.Tjdot = zero0
else:
    mod.Cadot = DerivativeVar(mod.Ca, initialize=-3.58709135E+01)
    mod.Tdot = DerivativeVar(mod.T, initialize=5.19191848E+03)
    mod.Tjdot = DerivativeVar(mod.Tj, initialize=-9.70467399E+02)
#: These guys as well (steady).
mod.Ca_ic = Param(mod.ncstr, default=1.9193793974995963E-02, mutable=True)
mod.T_ic = Param(mod.ncstr, default=3.8400724261199036E+02, mutable=True)
mod.Tj_ic = Param(mod.ncstr, default=3.7127352272578315E+02, mutable=True)

# m.Ca_ic = Param(m.ncstr, default=1.9193793974995963E-02)
# m.T_ic = Param(m.ncstr, default=3.8400724261199036E+02)
# m.Tj_ic = Param(m.ncstr, default=3.7127352272578315E+02)

mod.de_ca = Constraint(mod.t, mod.ncstr)
mod.de_T = Constraint(mod.t, mod.ncstr)
mod.de_Tj = Constraint(mod.t, mod.ncstr)

#: No need of these guys at steady.
if mod.steady:
    mod.Ca_icc = None
    mod.T_icc = None
    mod.Tj_icc = None
else:
    mod.Ca_icc = Constraint(mod.ncstr)
    mod.T_icc = Constraint(mod.ncstr)
    mod.Tj_icc = Constraint(mod.ncstr)

# let Ca0 := 1.9193793974995963E-02 ;
# let T0  := 3.8400724261199036E+02 ;
# let Tj0 := 3.7127352272578315E+02 ;

mod.kdef.rule = lambda m, i, n: _rule_k(m, i, n)
mod.de_ca.rule = lambda m, i, n: _rule_ca(m, i, n)
mod.de_T.rule = lambda m, i, n: _rule_t(m, i, n)
mod.de_Tj.rule = lambda m, i, n: _rule_tj(m, i, n)

if mod.steady:
    pass
else:
    mod.Ca_icc.rule = lambda m, n: _rule_ca0(m, n)
    mod.T_icc.rule = lambda m, n: _rule_t0(m, n)
    mod.Tj_icc.rule = lambda m, n: _rule_tj0(m, n)
    mod.Ca_icc.reconstruct()
    mod.T_icc.reconstruct()
    mod.Tj_icc.reconstruct()

mod.kdef.reconstruct()
mod.de_ca.reconstruct()
mod.de_T.reconstruct()
mod.de_Tj.reconstruct()

# Declare at framework level
# mod.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
# mod.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
# mod.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
# mod.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
# mod.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)
