#! env python
"""
Routines and interface for computing Mann-Whitney (Wilcoxon two-sample test)

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon

##
## Ported from perl code obtained from:
## http://www.fon.hum.uva.nl/Service/Statistics/Wilcoxon_Test.html

## COMMENTS are from original PERL code, which had the following
## information included in it:
##
##     Copyright (C) 1996, 2001  Rob van Son
##     
##     This program is free software; you can redistribute it and/or
##     modify it under the terms of the GNU General Public License
##     as published by the Free Software Foundation; either version 2
##     of the License, or (at your option) any later version.
##     
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##     
##     You should have received a copy of the GNU General Public License
##     along with this program(*); if not, write to the Free Software
##     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#############################################################################
## (*) Note: TAMO does not redistribute any source code from Rob van Son ####
#############################################################################

"""
from TAMO.util import Arith
from math import *
import sys

## Calculate k out n
def k_out_n(k,n):
    kn = 1.0
    while (k>0):
        kn *= float(n)/float(k)
        n = n-1
        k = k-1
    return kn

## This routine recursively counts the number of distributions of ranks over two
## samples for which the sum of the ranks in the smaller sample is smaller than or
## equal to a given upper bound W.
## W = the bound, Sum = the sum of ranks upto now, m-1 = one less than the
## number of elements in the smaller sample that still have to be done, 
## Start = the current position in the ranks list, *ankList = the array
## with all the ranks (this is NOT just the numbers from 1 - N because of ties).
## The list with ranks MUST be sorted in INCREASING order.
def CountSmallerRanks(W, Sum, m, Start, RankList):
    i, Temp, Smaller, End, mminus1 = 0, 0, 0, 0, 0
    if (Sum > W): return 0
    End = len(RankList)
    if (m > 0):
        mminus1 = m-1
        for i in range(Start,End-m):
            Temp = Sum + Ranklist[i]
            if (Temp > W): return Smaller
            Smaller += CountSmallerRanks(W, Temp, mminus1, i+1, Ranklist)
    else:
        if (Sum + End + 1 <= W): return (End - Start + 1);
        for i in range(Start,End):
            Temp = Sum + RankList[i]
            if (Temp <= W):
                Smaller += 1;
            else:
                return Smaller
    return Smaller

def main():
    ## Get input
    A = [float(x.strip()) for x in open(sys.argv[1]).readlines()]
    B = [float(x.strip()) for x in open(sys.argv[2]).readlines()]
    p, W = wlcxtest(A,B)
    print '%5.4g  ( %f )'%(p,W)

def WMWtest(A,B):
    """
    WMWtest(A,B) -- Computes the Wilcoxon-Mann-Whitney nonparametric W statistic for two distributions

    input:  list of numbers, list of numbers
    output: p-value, W-statistic
    """
    A.sort()
    B.sort()
    TotalList = A + B
    TotalList.sort()

    nA     = len(A)
    nB     = len(B)
    N      = nA + nB
    MaxSum = N*(N+1)/2.0
    H0     = MaxSum / 2.0
    
    
    ## Replace values by ranks
    previous = []
    start = 0
    Total_rank = TotalList[:]
    for i in range(len(TotalList)):
        if (TotalList[i] == previous):
            mean_rank = (start+i+2)/2.0
            for j in range(start,i+1):
                Total_rank[j] = mean_rank
        else:
            Total_rank[i] = i+1
            previous      = TotalList[i]
            start         = i
    
    ## Determine the shortest list
    if nA < nB: shortest = A
    else:       shortest = B
    nShortest = len(shortest);

    ## Summ the ranks in the shortest list
    W = 0
    for Value in shortest:
        i = 0
        while (i < len(TotalList) and Value != TotalList[i]): i += 1
        W += Total_rank[i]

    ## Use the smallest value of $W
    if (W > H0): W = MaxSum - W

    ## Determine the two-tailed level of significance
    p = 0

    ## First calculate the Normal approximation. This can be used to
    ## check whether a significant result is plausable for larger N.
    Permutations = k_out_n(nA, N)
    if (Permutations >= 25000) or (nShortest > 10):
        if W >= H0: Continuity = -0.5
        else:       Continuity =  0.5
        Z = (W+Continuity-nShortest*(N+1.0)/2.0)/sqrt(nA*nB*(N+1)/12.0);
        Z = fabs(Z)
        p = 2*(1-Arith.lzprob(Z))
    

    ## The exact level of significance, for large N, first check whether a
    ## significant result is plausable, i.e., the Normal Approximation gives
    ## a $p < 0.25.
    if (nShortest+1 < 10) and (p < 0.25) and (Permutations < 60000):
        # Remember that $W must be SMALLER than $MaxSum/2=$H0
        Less = CountSmallerRanks(W, 0 , len(shortest)-1, 0, Total_rank)
        # If $Less < $Permutations/2, we have obviously calculated the 
        # wrong way. We should have calculated UPWARD (higher than W)
        # We can't do that, but we can calculate $Less for $W-1 and
        # subtract it from $Permutations
        if (2*Less > Permutations):
            Less = CountSmallerRanks(W-1, 0, len(shortest)-1, 0, Total_rank)
            Less = Permutations - Less
        SumFrequencies = Permutations
        p = 2.0 * Less / SumFrequencies

    return p, W

if __name__ == '__main__': main()
