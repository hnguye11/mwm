from __future__ import division
import sys
from scipy.stats import expon, bernoulli
import matplotlib.pyplot as plt
from copy import copy
from math import exp, log10, log
import numpy as np
from time import time
from random import random, randint

sys.path.append('..')
from util import *


''' Implementation of Maximum Weight Minimization algorithm for Importance Sampling.
Example: Hypothetical Computer System Example (HCSE), "Fault tree handbook with aerospace applications", Stamatelatos 2002. '''


def CALC_LOGP(x, lbd, T):
    return sum([log(lbdi) - lbdi*xi if xi < T else -lbdi * T for xi,lbdi in zip(x,lbd)])


def CALC_LOGW(x, lbd, lbd1, T):
    return CALC_LOGP(x, lbd, T) - CALC_LOGP(x, lbd1, T)


def calc_f(x):
    HW,M1,M2,M3,M4,M5,MI1,MI2,P1,P2,PS,B1,B2 = x
    
    m1 = min(M1,MI1)
    m2 = min(M2,MI1)
    m3 = min(M3,max(MI1,MI2))
    m4 = min(M4,MI2)
    m5 = min(M5,MI2)
    MEM35 = sorted([m1,m2,m3,m4,m5])[2]

    IF = HW                     # suppose SW == OP == 0
    
    BUS = max([B1,B2])
    
    PROC1, PROC2 = P1, P2
    if PROC1 <= PROC2: PROC1 += PS
    else: PROC2 += PS        
    PROC = max(PROC1, PROC2)

    HECS = min([MEM35, IF, BUS, PROC])
    return HECS


def CALC_F(x):
    X = [x[int(i*M/N1):int((i+1)*M/N1)] for i in range(N1)]
    F = [IND(calc_f(x) >= T) for x in X]
    return IND(sum(F) < k)

    
def BINARY_SEARCH(x, y):
    assert(all([xi <= yi for xi, yi in zip(x,y)]))
    assert(CALC_F(x) == 1)
    
    if CALC_F(y) == 1: return y
    
    else:
        while max([yi-xi for xi,yi in zip(x,y)]) > 1:
            z = [int((xi+yi)/2) for xi,yi in zip(x,y)]
            if CALC_F(z) == 1: x = z
            else: y = z

        return x


def SAMPLE_RE_SET():
    while True:
        x = [randint(0,T) for i in range(M)]
        if CALC_F(x) == 1:
            break
        
    idx = np.random.permutation(M)
    
    for i in idx:
        y = copy(x)
        y[i] = T
        x = BINARY_SEARCH(x, y)

    return x


def APPX_RE_SET(filename, Ts=500, Ns=500):    
    X, Y = [], []
    start_time = curr_time = time()

    while curr_time - start_time < Ts and len(X) < Ns:
        x = SAMPLE_RE_SET()
        y = X_TO_GX(x)

        if not y in Y:
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


def X_TO_GX(x):
    y = []
    for i in range(M1):
        y += sorted([int(x[j]) for j in range(M) if EG_IDX[j]==i])
    return y


def GX_TO_X(gx):
    assert(len(gx) == M1)
    return [gx[i] for i in EG_IDX]


def NAIVE_MONTE_CARLO(Ns):
    F = [0] * Ns
    X = [expon(scale=1/lbdi).rvs(Ns) for lbdi in lbd]

    for ns in range(Ns):
        if ns % 10**4 == 0: print(ns)
        x = [Xi[ns] for Xi in X]
        F[ns] = CALC_F(x)
    
    
def IMPORTANCE_SAMPLING(q, Ns):
    lbd1 = [lbdi + qi for qi,lbdi in zip(q,lbd)]
    X = [expon(scale=1/lbdi).rvs(Ns) for lbdi in lbd1]
    W = [0] * Ns
    C = [0] * Ns
    
    for ns in range(Ns):
        if ns != 0 and ns % 10**4 == 0: print(ns)

        x = [Xi[ns] for Xi in X]
        W[ns] = exp(CALC_LOGW(x, lbd, lbd1, T))
        C[ns] = CALC_F(x)
        # print ns, W[ns], C[ns], x
        
    return W, C

##################################################

from hcse import *

##################################################

# APPX_RE_SET("Y_HCSE4_%d.csv"%k)
# exit()

##################################################

''' MWM solutions. '''
gq_4 = [0.001266, 0.000119, 0.000119, 0.000393, 0.000196, 0.000196, 0.000004]
gq_3 = [0.003170, 0.000218, 0.000218, 0.000646, 0.000297, 0.000297, 0.000007]
gq_2 = [0.005225, 0.000297, 0.000296, 0.000841, 0.000367, 0.000367, 0.000009]
gq_1 = [0.007373, 0.000363, 0.000363, 0.001009, 0.000423, 0.000423, 0.000011]

q = GX_TO_X(gq_2)               # parameter of the proposal distribution

Ns = 10**4
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
