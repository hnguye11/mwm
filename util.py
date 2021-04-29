from __future__ import division
from random import random as rnd, randint as rndi
import numpy as np
from math import sqrt


def IND(condition):
    return 1 if condition else 0

    
def READ_CSV_FILE(filename):
    data = []
    openfile = open(filename, "r")
    for line in openfile.readlines():
        tmp = line.strip().split(",")
        data += [[float(tmpi) for tmpi in tmp]]

    print("%s: %d entries"%(filename, len(data)))
    return data


def RANDOM(lb, ub):
    return lb + (ub - lb) * rnd()


def PROD(x):
    y = 1.0
    for xi in x: y *= xi
    return y


def TO_STRING(x, delim=", ", format="%.6f"):
    return delim.join([format%xi for xi in x])


def CALC_RELATIVE_ERROR(w):
    mean = np.mean(w)
    variance = np.var(w)
    std = sqrt(variance)
    return std / mean


''' Compute the effective sample size. '''
def CALC_ESS(w):
    if sum(w) == 0:
        return 0
    else:
        return sum(w)**2 / sum([wi**2 for wi in w])


def AVERAGE_EVERY(arr, k):
    return [np.mean(arr[i*k:(i+1)*k]) for i in range(int(len(arr)/k))]

