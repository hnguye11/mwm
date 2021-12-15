from __future__ import division
import sys
import networkx as nx
import numpy as np
from scipy.stats import bernoulli
from random import random, seed
from math import log10, sqrt
import matplotlib.pyplot as plt
from time import time

sys.path.append('..')
from util import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling. '''

def CALC_LOGP(x, q):
    return sum([log10(qi) if xi==1 else log10(1-qi) for xi,qi in zip(x,q)])


def CALC_LOGW(x, p, q):
    return CALC_LOGP(x, p) - CALC_LOGP(x, q)


def CALC_LOSS(x):
    G = nx.Graph()
    G.add_nodes_from(V)

    for i,Ei in enumerate(E):
        if x[i] == 1: G.add_edge(Ei[0], Ei[1])

    ''' Calculate S-T unreliability. '''
    rare_event = not nx.has_path(G, "s", "t")

    ''' RULE: return 1 if x is a rare event else 0. '''
    return 1 if rare_event else 0


def SAMPLE(q):
    ''' @param q: a proposal distribution of independent Bernoulli's. '''

    x = [1 if random() <= q[i] else 0 for i in range(M)]
    logw = CALC_LOGW(x, p, q)
    w = 10**logw
    c = CALC_LOSS(x)
    
    return x, w, c
    

''' Get a boundary sample from the rare event set (similar to generating D Spectra) '''
def SAMPLE_RE_SET():
    idx = np.random.permutation(range(M))
    x = [0] * M                 # start from a rare event sample

    for i in idx:
        x[i] = 1                # attempt to enable E[i] (for unreliability)
        if CALC_LOSS(x) == 0: x[i] = 0
        
    return x


''' Importance sampling using sampling vector q. '''
def IMPORTANCE_SAMPLING(q, Ns, plot=False):
    W = [0] * Ns
    C = [0] * Ns
    
    for ns in range(Ns):
        if ns != 0 and ns % 10**4 == 0: print(ns)
        x, W[ns], C[ns] = SAMPLE(q)

    if plot:
        logW = [log10(w) for w,c in zip(W,C) if c==1]
        print(len(logW), Ns)
        plt.hist(logW, bins=50, density=True, histtype="step")
        plt.show()
            
    return W, C
    

def X_TO_GX(x):
    gx = [0] * M1
    for i in range(M): gx[EG_IDX[i]] += x[i]
    return gx


def GX_TO_X(gx):    
    assert(len(gx) == M1)
    return [gx[i] for i in EG_IDX]


def NAIVE_MONTE_CARLO(Ns):
    W = [0] * Ns

    for ns in range(Ns):
        if ns % 1000 == 0: print(ns)
        x = [1 if random() < pi else 0 for pi in p]
        W[ns] = CALC_LOSS(x)

    return sum(W) / Ns


''' Sample the approximation of the rare event set using the order-
and the structure-preserving heuristic. 
@param Ts: time budget in seconds
@param Ns: sample budget.'''
def APPX_RE_SET(filename, Ts=500, Ns=500):
    X, Y = [], []
    start_time = curr_time = time()
    
    while curr_time - start_time < Ts and len(X) < Ns:
        x = SAMPLE_RE_SET()    # for unreliability
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
            
##################################################

from dodecahedron import *
# from lattice20x20 import *
M = len(E)
N = len(V)
Ns = 10**4
M1 = len(EG_UNIQUE)
p = p_3
gq = gq_3
q = GX_TO_X(gq)
assert(V[0] == "s")
assert(V[-1] == "t")
assert(len(p) == M)
SIM_INFO = "DODEC_M=%d_N=%d_Ns=%d"%(M, N, Ns)

print(SIM_INFO)
print("p  =", p)

##################################################

APPX_RE_SET("Y_dodec.csv", Ts=500, Ns=500)
exit()

##################################################

filename = "sim_result.txt"

''' Simulation. '''
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

