#!/usr/bin/env python3.6

from __future__ import division
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
from math import log10, exp
import matplotlib.pyplot as plt
import numpy as np
from time import time

from util import *
from util_lattice import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling.
Example: lattice networks in "Reliability Estimation for Networks withMinimum Flow Demand and Random Link Capacities", Botev 2018. '''


def UNPACK(u):
    q, t = UNZIP(u[:M1]), u[M1]
    return q, t


def UNZIP(gx):
    assert(len(gx) == M1)
    return [gx[i] for i in EG_IDX]

    
def f(u):
    global p, Logp
    
    q, t = UNPACK(u)
    Logq = [LOG(EXP_TWIST(p, q[i])) for i in range(M)]
        
    return [CALC_LOGW(x, Logp, Logq) - t for x in X]


def obj(u):
    return u[-1]


def solve(u0):
    # Setup
    u_bnds = [(-10, 0) for i in range(M1)] + [(-100,0)]
    f_lb = [-50] * len(X)
    f_ub = [10] * len(X)

    con = NonlinearConstraint(f, f_lb, f_ub)

    # Solve
    sol = minimize(obj, u0, method='SLSQP', bounds=u_bnds, constraints=con)

    return sol.x

##################################################

from lattice6x6 import *
# from lattice4x4 import *
eps = eps_8

X = READ_CSV_FILE("X_lattice6x6.csv")[:500]
X = [list(map(int, Xi)) for Xi in X]
print("len X = %d"%len(X))

p = GET_PROB(eps, rho, x_max)
Logp = [LOG(p)] * M

# initial guesses
u0 = [0] * M1 + [0]
start_time = time()
u1 = solve(u0)
end_time = time()

q, t = UNPACK(u1)
fu = f(u1)
print("q = [%s]"%TO_STRING(q))
print("t = %.3f"%t)
print("log(w) in [%f, %f]"%(min(fu) + t, max(fu) + t))
print("Elap = %.1f"%(end_time - start_time))
