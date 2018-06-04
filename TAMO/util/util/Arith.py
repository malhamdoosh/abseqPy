#Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
#All Rights Reserved
#
#Author: David Benjamin Gordon

import sys,math

def fact(a):
    '''
    fact(a)
    '''
    ans = 1.0
    while (a > 0):
        ans = ans * a
        a   = a   - 1
    return(ans)

        
def nlog10(x,min=1e-299):
    if x < min: x=min
    return math.fabs(math.log(x)/math.log(10))

def stirling(n):
    '''
    stirling(n)  (approximationg of log(n!))
    '''
    n = float(n)
    ans = (n + 0.5) * math.log(n) - n + 0.5*math.log(2*math.pi)
    return(ans)

def stircomb (x,y):
    '''
    stircomb (x,y)
    '''
    if (x < 100):    fx  = math.log(fact(x))
    else:            fx  = stirling(x)
    if (y < 100):    fy  = math.log(fact(y))
    else:            fy  = stirling(y)
    if (x-y < 100):  fxy = math.log(fact(x-y))
    else:            fxy = stirling(x-y)
    return(fx - fy - fxy)

def stirhypgeom (numinteresting, total, numpicked, numfound):
    '''
    stirhypgeom (numinteresting, total, numpicked, numfound)
    '''
    ans =   stircomb(numinteresting,         numfound) \
          + stircomb(total - numinteresting, numpicked - numfound) \
          - stircomb(total,                  numpicked)
    return(math.exp(ans))

def binomial(frac_exp,draws,hits):
    '''(frac,draws,hits) Binomial theorem expectation of getting exact number of hits'''
    frac_miss = 1-frac_exp
    pre_ans = stircomb(draws,hits)      +   \
              hits * math.log(frac_exp) +   \
              (draws-hits) * math.log(frac_miss)
    ans     = math.exp(pre_ans)
    #ans = math.exp(stircomb(draws,hits)) * \
    #      math.pow(frac_exp,  hits)      * \
    #      math.pow(frac_miss, draws-hits)
    #print math.exp(stircomb(draws,hits)), '*', \
    #      math.pow(frac_exp,  hits)      ,'*', \
    #      math.pow(frac_miss, draws-hits), '=', ans
    return ans

def binomialsumtail(frac_exp,draws,hits):
    '''(frac,draws,hits) Sum of the tail of the binomial distribution'''
    tot = 0.0
    if hits == 0: return 1.0
    epsilon = 1e-20
    for h in range(hits,draws+1):
        #print h, draws, stircomb(draws,h)
        try:
            newtot = tot + binomial(frac_exp,draws,h)
            if (10*h > draws) and ((newtot-tot)/(tot+epsilon) < 1e-5):  break
        except:
            sys.stderr.write('Problem evaluating %dth iteration of ( %d %d  ) [ %e ]   tot = %e '%(
                h,hits,draws,frac_exp,tot))
        tot = newtot
        #print tot
    return tot
    

def hypgeomsummore (numinteresting, total, numpicked, numfound):
    '''
    hypgeomsummore (numinteresting, total, numpicked, numfound)
    '''
    _min = min(numinteresting, numpicked)
    sum  = 0
    for i in range(numfound,_min+1):
        sum = sum + stirhypgeom(numinteresting,total,numpicked,i)
    return(sum)
    
def median(vals):
    N = float(len(vals))
    _v = vals[:]
    _v.sort()
    if N==1: return vals[0]
    i = int(N/2)
    if N/2 == i: #Average 
        median = 0.5 * (float(_v[i-1]) + float(_v[i]))
    else:
        median = float(_v[i])
    return median

def avestd(vals):
    (sum, sum2) = (0.,0.)
    N = float(len(vals))
    for val in vals:
        sum  = sum  + float(val)
        sum2 = sum2 + float(val)*float(val)
    if N == 0:
        ave = 0
        std = 0
    elif N == 1:
        ave = sum
        std = 0
    else:
        ave = sum /  N
        diff2 = sum2-(N*ave*ave)
        if diff2 >= 0.0:
            std = math.sqrt( diff2 / (N-1.0) )
        else:
            std = 0
    return(ave,std)


def norm_pvalue(ave,std,obs):
    return 1-lzprob((obs-ave)/std)

def lzprob(z):
    """
Returns the area under the normal curve 'to the left of' the given z value.
Thus, 
    for z<0, zprob(z) = 1-tail probability
    for z>0, 1.0-zprob(z) = 1-tail probability
    for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
Adapted from z.c in Gary Perlman's |Stat by Gary Strangman in stats.py
General Public License (GPL) v2, http://www.gnu.org/copyleft/gpl.html. 

Usage:   lzprob(z)
"""
    Z_MAX = 6.0    # maximum meaningful z-value
    if z == 0.0:
        x = 0.0
    else:
        y = 0.5 * math.fabs(z)
        if y >= (Z_MAX*0.5):
            x = 1.0
        elif (y < 1.0):
            w = y*y
            x = ((((((((0.000124818987 * w
                        -0.001075204047) * w +0.005198775019) * w
                      -0.019198292004) * w +0.059054035642) * w
                    -0.151968751364) * w +0.319152932694) * w
                  -0.531923007300) * w +0.797884560593) * y * 2.0
        else:
            y = y - 2.0
            x = (((((((((((((-0.000045255659 * y
                             +0.000152529290) * y -0.000019538132) * y
                           -0.000676904986) * y +0.001390604284) * y
                         -0.000794620820) * y -0.002034254874) * y
                       +0.006549791214) * y -0.010557625006) * y
                     +0.011630447319) * y -0.009279453341) * y
                   +0.005353579108) * y -0.002141268741) * y
                 +0.000535310849) * y +0.999936657524
    if z > 0.0:
        prob = ((x+1.0)*0.5)
    else:
        prob = ((1.0-x)*0.5)
    return prob

def rank_pvalue(obs,values):
    values.sort()   # 2 in [1, 2, 4, 7, 9, 12]
    if obs<0:
        values.reverse()
        obs = -obs   # 2 in [12, 9, 7, 4, 2, 1]
        remain = [x for x in values if (x < obs)]
    else:
        remain = [x for x in values if (x > obs)]

        
    numremain = float(max(1,len(remain)))

    p = numremain / len(values)

    #print numremain,len(values),p

    return p
    
