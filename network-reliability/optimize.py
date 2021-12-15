from __future__ import division
import sys
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
from math import log10
import matplotlib.pyplot as plt
from time import time

sys.path.append('..')
from util import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling. '''

def CALC_LOGP(x, p):
    return sum([log10(pi) if xi==1 else log10(1-pi) for xi,pi in zip(x,p)])


def CALC_LOGW(x, p, q):
    return CALC_LOGP(x, p) - CALC_LOGP(x, q)


def EXPAND(gx):
    assert(len(gx) == M1)
    return [gx[i] for i in EG_IDX]
    

def COMPRESS(x):
    assert(len(x) == M)
    gx = [0] * M1
    for i in range(M): gx[EG_IDX[i]] += x[i]
    return gx


def f(u):
    gq = u[:M1]
    q = EXPAND(gq)
    t = u[-1]
    return [CALC_LOGW(x, p, q) - t for x in X] + [sum(q) - sum(p)]
    

def obj(u):
    return u[-1]


def print_info(u):
    gq = u[:M1]
    q = EXPAND(gq)
    logw = [CALC_LOGW(x,p,q) for x in X]
    logp = [CALC_LOGP(x,q) for x in X]

    print("gq = [%s]"%TO_STRING(gq, format="%.6e"))
    print(min(logw), max(logw))
    print(sum(q) - sum(p))

##################################################

from dodecahedron import *

p = p_3

RELIABILITY = False
M = len(E)
N = len(V)
M1 = len(EG_UNIQUE)
gp = [p[EG_IDX.index(i)] for i in range(M1)]
X = READ_CSV_FILE("Y_dodec.csv")
assert(len(X[0]) == M)

if RELIABILITY:                 # increase edges probabilities
    BUDGET_MIN = 0
    BUDGET_MAX = M - sum(p)
else:                           # decrease edge probabilities
    BUDGET_MIN = -sum(p)
    BUDGET_MAX = 0

##################################################

f_lb = [-1000] * len(X) + [BUDGET_MIN]
f_ub = [50] * len(X) + [BUDGET_MAX]

if RELIABILITY:
    u_bnds = [(pi, 1-1e-6) for pi in gp] + [(-100,0)]
else:
    u_bnds = [(1e-6, pi) for pi in gp] + [(-100,0)]

start_time = time()

while True:
    ''' Initial guess '''
    # u0 = [pi+0.01 for pi in gp] + [0]
    if RELIABILITY:
        u0 = [RANDOM(pi, 1) for pi in gp] + [0]
    else:
        u0 = [RANDOM(1e-6, pi) for pi in gp] + [0]
    f0 = f(u0)

    if all([f1<=f2<=f3 for f1,f2,f3 in zip(f_lb,f0,f_ub)]):
        ''' Boundary conditions satisfied! '''
        print("\nInitial guess:")
        print_info(u0)
        break

con = NonlinearConstraint(f, f_lb, f_ub)
sol = minimize(obj, u0, method='SLSQP', bounds=u_bnds, constraints=con)
u1 = sol.x

end_time = time()
elap_time = end_time - start_time
print("Solution:")
print_info(u1)
print("Elap = %.1f"%elap_time)

