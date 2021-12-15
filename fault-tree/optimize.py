from __future__ import division
import sys
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
import matplotlib.pyplot as plt
import numpy as np
from time import time
from math import log, exp

sys.path.append('..')
from util import READ_CSV_FILE, TO_STRING

''' Fault tree handbook with aerospace applications. 
The Hypothetical Computer System Example (HCSE). '''


def CALC_LOGP(x, lbd, T):
    return sum([log(lbdi) - lbdi*xi if xi < T else -lbdi * T for xi,lbdi in zip(x,lbd)])


def CALC_LOGW(x, lbd, lbd1, T):
    return CALC_LOGP(x, lbd, T) - CALC_LOGP(x, lbd1, T)


def UNPACK(u):
    q, t = GX_TO_X(u[:M1]), u[M1]
    return q, t



def GX_TO_X(gx):
    assert(len(gx) == M1)
    return [gx[i] for i in EG_IDX]

    
def f(u):
    global lbd, T
    
    q, t = UNPACK(u)
    lbd1 = [lbdi + qi for qi,lbdi in zip(q,lbd)]
        
    return [CALC_LOGW(x, lbd, lbd1, T) - t for x in X]


def obj(u):
    return u[-1]


def solve(u0):
    # Setup
    u_bnds = [(0, 1) for i in range(M1)] + [(-100,0)]
    f_lb = [-50] * len(X)
    f_ub = [10] * len(X)

    con = NonlinearConstraint(f, f_lb, f_ub)

    # Solve
    sol = minimize(obj, u0, method='SLSQP', bounds=u_bnds, constraints=con)

    return sol.x

##################################################

lbd_hw = 5e-5                   # Hardware
lbd_mi = 5e-5                   # Memory Interface
lbd_m = 6e-5                    # Memory
lbd_p = 1e-4                    # Processor
lbd_ps = 1e-4                   # Processor Swap
lbd_b = 1e-6                    # Bus
N1 = 4                          # Number of replicas
k = 2                           # Mininum number of working replicas
T = gmm = 100                   # Mission time
lbd = [lbd_hw] + [lbd_m] * 5 + [lbd_mi] * 2 + [lbd_p] * 2 + [lbd_ps] + [lbd_b] * 2
EG_IDX = [0, 1, 1, 2, 1, 1, 3, 3, 4, 4, 5, 6, 6] 
lbd *= N1
EG_IDX *= N1
M = len(lbd)
M1 = len(set(EG_IDX))

##################################################

X = READ_CSV_FILE("Y_HECS4_%d.csv"%k)[:500]
X = [list(map(int, Xi)) for Xi in X]
print("len X = %d"%len(X))

# initial guesses
u0 = [0] * M1 + [0]
start_time = time()
u1 = solve(u0)
end_time = time()
q, t = UNPACK(u1)
fu = f(u1)
print("q_zip = [%s]"%TO_STRING(u1[:M1]))
print("t = %.3f"%t)
print("log(w) in [%f, %f]"%(min(fu) + t, max(fu) + t))
print("elap = %.1f"%(end_time-start_time))
