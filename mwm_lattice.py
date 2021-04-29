from __future__ import division
import sys
import networkx as nx
from copy import copy
import numpy as np
from scipy.stats import rv_discrete
from math import log10, sqrt
import matplotlib.pyplot as plt
from time import time

from util import *
from util_lattice import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling.
Example: lattice networks in "Reliability Estimation for Networks withMinimum Flow Demand and Random Link Capacities", Botev 2018. '''


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


def SAMPLE_RE():
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


def APPX_RARE_EVENT_SET(filename, Ts=500, Ns=500):
    X, Y = [], []
    start_time = curr_time = time()

    while curr_time - start_time < Ts and len(X) < Ns:
        x = SAMPLE_RE()
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
            print len(Y), y, "%.1f, %.1f"%(elap, total_elap)
    
    
def IMPORTANCE_SAMPLING(q, Ns):
    p = GET_PROB(eps, rho, x_max)
    Logp = [LOG(p)] * M
    Q = [EXP_TWIST(p, q[i]) for i in range(M)]
    Logq = [LOG(Qi) for Qi in Q]        
    Rv = [rv_discrete(values=(range(x_max+1),Qi)).rvs(size=Ns) for Qi in Q]
    W = [0] * Ns
    C = [0] * Ns
    
    for ns in range(Ns):
        if ns != 0 and ns % 10**3 == 0: print ns

        x = [Rvi[ns] for Rvi in Rv]
        W[ns] = 10**CALC_LOGW(x, Logp, Logq)
        C[ns] = IND(CALC_F(x) < gmm)
        
    return W, C

##################################################

from lattice4x4 import * 
# from lattice6x6 import *
eps = eps_8
q = q_8

##################################################

APPX_RARE_EVENT_SET("X_lattice6x6.csv")
exit()

##################################################

Ns = 10**4
filename = "mu_lattice4x4_8.txt"

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

    print "est : ", mu_is
    print "hit : ", hit, "(%.2f%%)"%(hit * 100 / Ns)
    print "max : ", mu_max
    print "re  : ", re
    print "ess : ", ess
    print "elap: ", elap
    print "---------------"

    openfile = open(filename, "a")
    openfile.write("%d, %.8e, %d, %.5f, %.5f, %.5f\n"%(Ns, mu_is, hit, mu_max, re, ess))
    openfile.close()

exit()
