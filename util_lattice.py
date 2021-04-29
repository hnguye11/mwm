from __future__ import division
from math import log10, exp


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
