from __future__ import division
import sys
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
from math import log10, exp
import matplotlib.pyplot as plt
import numpy as np
from time import time

sys.path.append('..')
from util import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling.
Example: lattice networks in "Reliability Estimation for Networks with Minimum Flow Demand and Random Link Capacities", Botev 2018. '''


def LOG(p):
    return [log10(pi) for pi in p]


def GET_PROB(eps, rho, x_max):
    p = [eps * rho**(x_max-xi-1) for xi in range(x_max)]
    p += [1 - sum(p)]
    return p


def CALC_LOGP(x, Logp):
    assert(len(x) == len(Logp))
    return sum([logpi[xi] for xi,logpi in zip(x,Logp)])


def CALC_LOGW(x, Logp, Logq):
    return CALC_LOGP(x, Logp) - CALC_LOGP(x, Logq)


def EXP_TWIST(p, theta):
    q = [exp(i*theta) * pi for i,pi in enumerate(p)]
    q_sum = sum(q)
    return [qi/q_sum for qi in q]


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

# from lattice6x6 import *
from lattice4x4 import *
eps = eps_8

X = READ_CSV_FILE("Y_lattice4x4.csv")
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
