"""
Poisson.py

Utilities for estimating and computing Poisson distributions.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""
import sys
from math  import *
from Arith import fact, avestd, stirling

def Poisson_estimate(obs,_dist):
    '''
    Poisson_estimate(obs,_dist) -- Compute optimal lambda for input distribution
    and return P and the distance between the fit curve and the input distribution.
    '''
    bestlam,bestdist = bestPoissondist(_dist)
    #print bestlam,bestdist
    #sys.stdout.flush()
    p = Poisson_sumtail(int(obs),bestlam)
    return p, bestdist


def Poisson(k,lam):
    '''
    Poisson(k,lam)
    '''
    try:
        return exp(-lam)*pow(lam,k)/fact(k)
    except:  pass

    try:
        lans = -lam + (k*log(lam)-stirling(k))
        return exp(lans)
    except:
        #return 0.5
        print "# Bad Values for Poisson(k,lam) =: ",k,lam
        raise "BadPoissonParams"
        #sys.stdout.flush()
        #sys.exit(1)

def Poisson_sumtail(k,lam):
    '''
    Poisson_sumtail(k,lam) -- Sum up the tail (P-value)
    '''
    total = 0.0
    i = k
    while 1:
        #print total,i,k,lam
        try:
            p = Poisson(i,lam)
        except BadPoissonParams:
            if total > 0.9:
                total = 1.0
                break
            else:
                print "# Bad Values for Poisson(i,lam) tot=: ",i,lam,total
                raise 
        if p < (total/1e10):
            break
        total += p
        i += 1
    return total

def bestPoissondist(_dist):
    '''
    bestPoissondist(_dist)  -- Normalize and find the best matching Poisson distribution
    '''
    D = {}
    total = float(len(_dist))

    for n in _dist:
        try:    D[int(n)] += 1
        except: D[int(n)] = 1
    maxval = max(D.keys())+2 #Include a "zero" in the tail
    xvals = [i for i in range(maxval)]
    yvals = []
    for i in range(maxval):
        if D.has_key(i): yvals.append(D[i]/total)
        else:            yvals.append(0)
    return bestPoisson(_dist,xvals,yvals)

def bestPoisson(_dist,xvals,yvals): 
    '''
    bestPoisson(_dist,xvals,yvals)  -- Given a curve described as a set of (x,y) pairs,
    find the value of lambda that provides the best fit Poisson distributon, and return
    the distance between the distributions as well.
    '''
    pairs    = zip(xvals,yvals)
    ave,std  = avestd(_dist)
    var      = std*std
    minlam   = min(ave,var)*0.9
    maxlam   = max(ave,var)*1.1
    interval = (maxlam-minlam)/100.0
    #print ave,std,var,minlam, maxlam
    sys.stdout.flush()

    dstats = []
    for i in range(100):
        lam = minlam + interval*i
        dtot = 0
        for x,y in pairs:
            pred  = Poisson(x,lam)
            dtot += (y-pred)*(y-pred)
        #RMS = sqrt(dtot/len(pairs))
        dstats.append((dtot,lam))

    dstats.sort()
    dbest, lambest = dstats[0]
    #RMSbest, lambest = dstats[0]

    ymax = float(max(yvals))
    scale = 60 / ymax

    old_fit = fabs(ave-var)/ave

    HISTOGRAM = 0
    if HISTOGRAM:
        print "##Over interval %f - %f: best lambda = %s  d= %f  old: %7.4f"%(
            minlam,maxlam,lambest,dbest,old_fit)

        for x,y in pairs:
            pred = Poisson(x,lambest)
            txtheight = int(scale*y+0.5)
            line = 'p'*txtheight + ' '*(80-txtheight)
            baridx = int(scale*pred+0.5)
            if baridx !=0:
                line = line[0:baridx] + '|' + line[baridx+1:]
            print '## %4d %7.5f %7.5f %s'%(x, y, pred, line)

    #return lambest,RMSbest
    return lambest,dbest
