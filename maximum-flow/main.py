from __future__ import division
import sys
import networkx as nx
import numpy as np
from scipy.stats import rv_discrete
from copy import copy
from math import log10, sqrt, exp
import matplotlib.pyplot as plt
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


def CALC_F(x):
    G = nx.DiGraph()
    G.add_nodes_from(V)

    for i,Ei in enumerate(E):
        G.add_edge(Ei[0], Ei[1], capacity=x[i])

    res = nx.maximum_flow(G, "s", "t")
    return res[0]


def BINARY_SEARCH(x, y):
    assert(all([xi <= yi for xi, yi in zip(x,y)]))
    assert(CALC_F(x) < gmm)
    
    if CALC_F(y) < gmm: return y
    
    else:
        while max([yi-xi for xi,yi in zip(x,y)]) > 1:
            z = [int((xi+yi)/2) for xi,yi in zip(x,y)]
            if CALC_F(z) < gmm: x = z
            else: y = z

        return x


def SAMPLE_RE_SET():
    while True:
        x = [rndi(0,x_max) for i in range(M)]
        if CALC_F(x) < gmm: break
        
    idx = np.random.permutation(M)
    
    for i in idx:
        y = copy(x)
        y[i] = x_max
        x = BINARY_SEARCH(x, y)

    return x


def X_TO_GX(x):
    y = []
    for i in range(M1):
        y += sorted([int(x[j]) for j in range(M) if EG_IDX[j]==i])
    return y


def APPX_RE_SET(filename, Ts=500, Ns=500):
    X, Y = [], []
    start_time = curr_time = time()

    while curr_time - start_time < Ts and len(X) < Ns:
        x = SAMPLE_RE_SET()
        y = X_TO_GX(x)

        if not x in X and not y in Y:
            X.append(x)
            Y.append(y)

            openfile = open(filename, "a")
            openfile.write(",".join(map(str, x)) + "\n")
            openfile.close()

            next_time = time()
            elap = next_time - curr_time
            total_elap = next_time - start_time
            curr_time = next_time
            print(len(Y), y, "%.1f, %.1f"%(elap, total_elap))
    
    
def IMPORTANCE_SAMPLING(q, Ns):
    p = GET_PROB(eps, rho, x_max)
    Logp = [LOG(p)] * M
    Q = [EXP_TWIST(p, q[i]) for i in range(M)]
    Logq = [LOG(Qi) for Qi in Q]        
    Rv = [rv_discrete(values=(range(x_max+1),Qi)).rvs(size=Ns) for Qi in Q]
    W = [0] * Ns
    C = [0] * Ns
    
    for ns in range(Ns):
        if ns != 0 and ns % 10**4 == 0: print(ns)

        x = [Rvi[ns] for Rvi in Rv]
        W[ns] = 10**CALC_LOGW(x, Logp, Logq)
        C[ns] = IND(CALC_F(x) < gmm)
        
    return W, C

##################################################

from lattice4x4 import * 
# from lattice6x6 import *
Ns = 10**4
eps = eps_8
p = GET_PROB(eps,rho,x_max)
q = q_8
assert(V[0] == "s")
assert(V[-1] == "t")
SIM_INFO = "LATTICE_M=%d_N=%d_Ns=%d"%(M, N, Ns)

print(SIM_INFO)
print("p  =", p)

##################################################

# APPX_RE_SET("Y_lattice4x4.csv")
# exit()

##################################################

filename = "sim_result.txt"

while True:
    start_time = time()
    W, C = IMPORTANCE_SAMPLING(q, Ns)
    elap = time() - start_time
    Mu = [w*c for w,c in zip(W,C)]    
    mu_is = sum(Mu) / Ns
    hit = np.count_nonzero(Mu)
    re = CALC_RELATIVE_ERROR(Mu)
    ess = CALC_ESS(Mu)
    mu_max = log10(max(Mu)) if hit > 0 else -100

    print("est : ", mu_is)
    print("hit : ", hit, "(%.2f%%)"%(hit * 100 / Ns))
    print("max : ", mu_max)
    print("re  : ", re)
    print("ess : ", ess)
    print("elap: ", elap)
    print("---------------")

    openfile = open(filename, "a")
    openfile.write("%d, %.8e, %d, %.5f, %.5f, %.5f\n"%(Ns, mu_is, hit, mu_max, re, ess))
    openfile.close()

exit()
