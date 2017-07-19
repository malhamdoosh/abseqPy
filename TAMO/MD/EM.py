'''
EM.py

This file contains class definitions for 5 classes useful for using EM for
Motif Discovery


class EM:                 Implementation of an EM

class MarkovBackground:   Implementation of a n-th order Markov Background model
                          (Pre-compute models using Background.py)

class Probe:              A probe (string with additional baggage)

class MotifCandidate:     Essentially a "super" Motif, that can compute its MAP score
                          probably should have been subclassed.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
'''

import sys, re, os, math, time 
# import Numeric   #Only if using MotifCandidate._plogp_new 
import time
import MDsupport

import TAMO.paths
from   TAMO    import MotifTools
from   TAMO.seq import Fasta


theMarkovBackground = None

class MotifCandidate:
    '''
    A candidate for Motif consists of the following:
       A set of segments "wmers"
       A pssm computed from these segments
       A score (probably the MAP score described in the literature)
       General propeties:
         # segments
         Width

       A candidate must also be able to:
         Evaluate its own score
         Modify itself:
           Add wmer
           Remove wmer      
    '''
    
    def __init__(self, wmers=''):
        self.MAP  = 0
        self.log2 = math.log(2)
        self.Q = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast defaults
        self.last_bgprob = 0
        self.deltaMAP = 0.0

        global theMarkovBackground
        if theMarkovBackground:
            self.bgprob = theMarkovBackground.zeroth()

        if wmers:
            self.wmers    = wmers
            self._update()
            self.needs_update = 1

    def find_wmers(self,seqs):
        self.wmers = []
        for seq in seqs:
            (bestscore,bestmatch) = (0,'')
            (matches,endpoints,scores) = self.pssm._scan(seq,-10)
            for match,score in zip(matches,scores):
                if score > bestscore:
                    bestscore = score
                    bestmatch = match
            if bestmatch:
                self.wmers.append(bestmatch)
        self._update()

    def _plogp_new(self,wmers): #NOT USED (Also, reuquires Numeric)
        nwmers = float(len(wmers))
        toAscii = {'A':97, 'C':99, 'G':103, 'T':116}
        sumsD = {}
        AsciiSeq = Numeric.array(wmers).astype(type(1))
        plogp = 0
        for key in toAscii.keys():
            sumsD[key] = Numeric.sum(Numeric.equal(AsciiSeq,toAscii[key])).astype(type(1.))
            f = sumsD[key]/nwmers + 0.0000000001
            plogp = plogp + Numeric.sum(Numeric.log(f) * f / Numeric.log(2))
        return(plogp)

    def _plogp(self,wmers):
        nwmers = float(len(wmers))
        plogp  = 0
        for i in range(len(wmers[0])):
            C = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N':0}
            for wmer in wmers:
                C[wmer[i]] = C[wmer[i]] + 1
            for key in C.keys():
                if key == 'N': continue
                f = C[key]/nwmers
                if f > 0:
                    plogp = plogp + f * math.log(f)/self.log2
        return(plogp)

    def _plogp_memo(self, wmer, count):
        nwmers = float(len(self.wmers) + count)
        plogp = 0
        for i in range(len(wmer)):
            letter = wmer[i]
            for key in ['A', 'C', 'T', 'G']:
                lettercount = self.pssm.counts[i][key]
                if key == letter:
                    lettercount = lettercount + count
                #print lettercount, nwmers
                f = float(lettercount)/nwmers
                if f > 0:
                    plogp = plogp + f * math.log(f)/self.log2
        return(plogp)


    def computeMAP_memo(self, wmer, count):
        wmers = self.wmers
        nwmers = float(len(wmers))  + count
        width  = float(len(wmers[0]))
        plogp  = self._plogp_memo(wmer,count)
        bgprob = self.last_bgprob + count * self.logPbackground(wmer)
        #print "MEMO: %d %f %f"%(nwmers,plogp,bgprob)
        MAP = math.log(nwmers)/ (self.log2 * width) \
              * (plogp - (1.0/nwmers) * bgprob)
        #MAP = nwmers/ width * (plogp - (1.0/nwmers) * bgprob - math.log(75)/math.log(2))
        return(MAP,bgprob)
        
    def _update(self):
        self.nwmers   = len(self.wmers)
        self.pssm     = MotifTools.Motif(self.wmers,self.Q)
        self.MAP      = self.computeMAP() 
        self.pssm.MAP = self.MAP
        
    def _recompute(self,_seqs):
        self.nwmers       = len(_seqs)
        self.wmers        = _seqs
        self.pssm         = MotifTools.Motif(_seqs,self.Q)
        #self.last_bgprob  = _bgprob
        self.MAP          = self.computeMAP()
        self.pssm.MAP     = self.MAP
        self.needs_update = 1

    def computeMAP(self,in_wmers=''):
        if in_wmers:
            wmers = in_wmers
        else:
            wmers = self.wmers
        nwmers = float(len(wmers))
        width  = float(len(wmers[0]))
        plogp  = self._plogp(wmers)
        bgprob = 0
        for wmer in wmers:
            bgprob = bgprob + self.logPbackground(wmer)
        #print "ORIG: %d %f %f"%(nwmers,plogp,bgprob)
        MAP = math.log(nwmers)/ (self.log2 * width) * (plogp - (1.0/nwmers) * bgprob)
        if not in_wmers:
            self.last_bgprob = bgprob
        return(MAP)
    
    def logPbackground(self,wmer):
        if theMarkovBackground:
            Ptot = theMarkovBackground.logbackground(wmer)
        else:
            Ptot = math.log(0.25)/self.log2 * len(wmer)
        return(Ptot)
    def __repr__(self):
        s = ''
        s = s + '%-10s  (Bits: %6.2f   MAP: %6.2f   D: %5.2f %2d)'%\
                 (self.pssm.__repr__(), self.pssm.totalbits,self.MAP,
                  self.pssm.seeddist,self.pssm.seednum)
        return s

    def check_and_update(self,wmer,count,verbose=''):
        if self.has_wmer(wmer):
            return
        scanscore,maxscore = self.pssm._scan(wmer,-1000,"FORW_ONLY")[2][0],self.pssm.maxscore
        if (scanscore/maxscore) < 0:  #This is a shortcut: Use Scan to decide
            return()                  #Whether the full MAP calculation is worthwhile
        
        # Does adding one copy help?  (we don't want to overwhelm with
        # all copies.)
        #_wmers = self.wmers[:]
        #_wmers.append(wmer)
        #_MAP  = self.computeMAP(_wmers)
        (_MAP,bgprob) = self.computeMAP_memo(wmer,1)
        #print ('%8.4f %8.4f %8.4f %8.4f ')%(
        #    self.MAP,_MAP,self.computeMAP(_wmers), \
        #    _MAP-self.computeMAP(_wmers))
        delta = _MAP - self.MAP
        if 0:  #Used for printing scanscore vs MAP score stats
            print '## %10.5f %10.5f %10.5f   %10.5f %10.5f %10.5f'%(
                scanscore,maxscore,scanscore/maxscore,_MAP,self.MAP,_MAP/self.MAP),
            if _MAP > self.MAP:
                print 1,wmer
            else:
                print 0,wmer
            return()
        
        if _MAP > self.MAP + self.deltaMAP:
            # The one wmer helps the map score.  So add all its copies
            # This seems like the right thing to do... doesn't it?  How
            # could we justify only adding _some_ of the copies?
            _wmers = self.wmers[:]
            for i in range(count):
                _wmers.append(wmer)
            _pssm = MotifTools.Motif(_wmers,self.Q)
            (_MAP,self.last_bgprob)  = self.computeMAP_memo(wmer,count)
            if verbose: print 'Adding %d copies of %s: O:%f N:%f (delta:%f)'% \
               (count,wmer,self.MAP,_MAP,delta)
            sys.stdout.flush()
            self.nwmers  = len(_wmers)
            self.wmers   = _wmers
            self.pssm    = _pssm
            self.MAP     = _MAP
            self.pssm.MAP= _MAP
            self.needs_update = 1
            return("ACCEPTED")

    def MAPscan(self,nmers):
        (maxMAP,maxnmer) = (0,'')
        for nmer in nmers:
            scanscore = self.pssm._scan(nmer,-1000,"Forw_Only")[2][0]
            if scanscore < 0: continue
            (_MAP, bgprob) = self.computeMAP_memo(nmer,1)
            if _MAP > maxMAP:
                maxMAP  = _MAP
                maxnmer = nmer
        return(maxnmer)

    def MAPpurge(self,verbose=''):
        Tot_Removed = 0
        while 1:
            wmers    = self.wmers
            new_list = self._MAPpurge_list()
            delta = len(wmers) - len(new_list)
            Tot_Removed = Tot_Removed + delta
            if delta == 0:  #Nothing purged
                break
            if verbose:
                print "\tPurged %4d (Total %4d) sequences (leaving %4d) %s"% \
                      (delta,Tot_Removed,len(new_list),self.pssm.oneletter)
                sys.stdout.flush()
            self._recompute(new_list)
            #print self

    def _MAPpurge_list(self):
        omit_list = []
        for wmer in self.wmers:
            count = 0
            for _wmer in self.wmers:
                if wmer == _wmer:
                    count = count + 1
            if count == len(self.wmers):  #Only one type of wmer in this thing
                return(self.wmers)        #Don't purge it!!!
            (_MAP,_bgprob) = self.computeMAP_memo(wmer,-count)
            if _MAP > self.MAP:
                omit_list.append(wmer)
        new_list  = []
        for wmer in self.wmers:
            if not (wmer in omit_list):
                new_list.append(wmer)
        return(new_list)

    def purge(self,verbose = ''):
        if verbose: print "Purging %s"%self
        #Build list of wmers
        wmersD = {}
        for wmer in self.wmers:
            if not wmersD.has_key(wmer):
                wmersD[wmer] = 1
        for wmer_omit in wmersD.keys():
            #Build temporary sequence list
            _seqs = []
            for wmer in self.wmers:
                if wmer != wmer_omit:
                    _seqs.append(wmer)
            count = len(self.wmers) - len(_seqs)
            if len(_seqs) == 0: continue
            (_MAP,_bgprob)  = self.computeMAP_memo(wmer_omit,-count)
            #print ('%8.4f %8.4f %8.4f %8.4f ')%(
            #    self.MAP,_MAP,self.computeMAP(_seqs), \
            #    _MAP-self.computeMAP(_seqs))
            #_MAP  = self.computeMAP(_seqs)
            if _MAP > self.MAP:
                if verbose: print 'Purging %d copies of %s: O:%f N:%f (delta:%f)'% \
                   (count,wmer_omit,self.MAP,_MAP,_MAP-self.MAP)
                sys.stdout.flush()
                self.nwmers  = len(_seqs)
                self.wmers   = _seqs
                self.pssm    = MotifTools.Motif(_seqs,self.Q)
                self.last_bgprob = _bgprob
                self.MAP     = _MAP
                self.pssm.MAP= _MAP
                self.needs_update = 1
        
    def has_wmer(self,wmer):
        rc = MotifTools.revcomplement(wmer)
        if (wmer in self.wmers) or (rc in self.wmers):
            return(1)
        else:
            return(0)

class EM:
    def __init__(self,seed_seqs, all_seqs, width = 6, verbose = ''):
        self.seed_seqs  = seed_seqs #Sequences to be scanned for seeds
        self.seqs       = all_seqs
        self.candidates = []
        self.models     = []      #Set directly or computed from seed_seqs
        self.width      = width
        self.verbose    = verbose
        if width:
            self.goodwmersT = MotifTools.top_nmers(self.width,self.seed_seqs,1,"")
        else:
            self.goodwmersT = zip(self.seed_seqs,range(len(self.seed_seqs)))
        self.bgprob     = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31}
        self.beta       = 0.001
        self.deltamin   = 1e-3
        self.probes     = []
        self.method     = "ZOOPS" # OOPS or ZOOPS )
        self.param      = {}
        self.gapflank   = 0
        self.gapweight  = 0.2
        self.seedbeta   = 0.02
        self.joint      = 1

        global theMarkovBackground
        if theMarkovBackground:
            self.bgprob = theMarkovBackground.zeroth()

        '''DELETE
        if all_seqs:  #Should we Go?
            self.EM()
        '''

    def report(self):
        if self.width:
            print "# A total of %5d %d-mer seeds specified or found in top %3d sequences"% \
                  (len(self.goodwmersT),self.width,len(self.seed_seqs))
        else:
            print "# A total of %5d seeds specified or found in top %3d sequences"% \
                  (len(self.goodwmersT),len(self.seed_seqs))
            
        print "# Using beta = %f   and  Converging once d(theta) < %f"% \
              (self.beta, self.deltamin)
        if theMarkovBackground:
            print "# Using High-order (%d-th)  Markov Background model"% \
                  (theMarkovBackground.highestorder -1)
        if self.gapflank > 0:
            print "# Using Gapped model: flank %d, weight %3.1f"%(self.gapflank, self.gapweight)
        sys.stdout.flush()

    def seed_models(self):
        V = self.verbose
        beta = self.seedbeta
        for nmer,count in self.goodwmersT:
            M = MotifTools.Motif(None,self.bgprob) #No automatic computation
            M.compute_from_nmer(nmer,beta)
            self.models.append(M)

    def calcmask(self,width):
        self.mask = map(lambda x:1,range(width))
        if self.gapflank > 0:
            gf = self.gapflank
            gw = self.gapweight
            for pos in range(gf,width-gf):
                self.mask[pos] = gw
            self.mask[0]  = gw
            self.mask[-1] = gw

    def EM_Cstart(self):
        verbose = self.verbose
        if verbose:
            print "Seeding models..."
            sys.stdout.flush()
        self.seed_models()

        #Initialize parameters
        if not self.param.has_key('gamma'): self.param['gamma'] = 0.2
        timings = {'Probes':0, 'Background':0, 'C EM':0, 'Post':0}
        _time = time.time()

        for seq in self.seqs:
            P = Probe(seq)
            self.probes.append(P)


        _time2 = time.time(); timings['Probes'] = _time2-_time; _time = _time2

        if verbose: print "Optimizing candidates by EM."
        if verbose: sys.stdout.flush()

        c_logZ_sets = {}
        for Model,i in zip(self.models,range(len(self.models))):
            width = Model.width

            self.calcmask(width)
              
            if not c_logZ_sets.has_key(width):
                c_logZs_set = []
                if verbose: print "#%s   |%s|"%(' '*28,'-'*len(self.seqs))
                if verbose: sys.stdout.flush()
                if verbose: print "Computing background (width %2d)  "%width,
                for P in self.probes:
                    if verbose: sys.stdout.write('.')
                    if verbose: sys.stdout.flush()
                    logZs = self.all_Wmers(width,P)
                    c_logZs = MDsupport.list2double(logZs)
                    c_logZs_set.append(c_logZs)
                    #P.c_wmerbgs = MDsupport.list2double(map(lambda x: x.logQtot, Wlist))
                c_logZ_sets[width] = c_logZs_set
                if verbose: print

            c_logZ_set = c_logZ_sets[width]
            for P,c_logZs in zip(self.probes,c_logZ_set):
                P.c_wmerbgs = c_logZs
                
            _time2 = time.time()
            timings['Background'] = timings['Background'] +_time2-_time
            _time = _time2


            '''Perform EM'''
            _time  = time.time()
            newModel = self.EM_C(Model, self.probes)
            _time2 = time.time(); timings['C EM'] = timings['C EM'] + _time2-_time; _time = _time2

            #print "cLL: ",newModel.joint
            #print "pLL: ",self.compute_joint(newModel,Wmers_by_seq)

            '''Was there a problem?'''
            if newModel == None:
                continue


            '''Set various things in PSSM'''
            #Distance(s)
            seeddist = MotifTools.infomaskdiff(newModel,Model)
            print '%s ----> %s'%(Model,newModel)
            print "Seed %2d: %s  -->  %s  mask:%9.5f  infoMask:%9.5f d:%9.5f"%(
                i, Model, newModel,
                MotifTools.maskdiff(newModel,Model),
                MotifTools.infomaskdiff(newModel,Model), #order is important
                Model-newModel)
            #Seed
            if Model.seedtxt: newModel.seedtxt = Model.seedtxt
            if Model.source:  newModel.source  = Model.source
            
            #newModel.denoise()
            newModel.seeddist = seeddist
            newModel.seednum  = i
            print newModel
            newModel._print_p()
            newModel._print_ll()

            '''Set various things in Candidate (like a wrapper for PSSM)'''
            C = MotifCandidate()
            C.pssm = newModel.copy()
            #C.wmers = self.best_by_Z(Wmers_by_seq)
            C.wmers  = [newModel.emit() for junk in range(20)]
            #C._update()  #MAJOR REMOVAL????????? DBG 10-14-03
            #C.MAPpurge()
            C.pssm = newModel.copy()  
            self.candidates.append(C)
            _time2 = time.time(); timings['Post'] = timings['Post']+_time2-_time;_time = _time2

        '''Print Timing Information'''
        if verbose:
            print "# Timing Information"
            _t = 0
            for timing in timings.keys():
                _t = _t + timings[timing]
            for timing in timings.keys():
                print "# %12s %f  %f%%"%(timing,timings[timing],timings[timing]*100/_t)


    def compute_joint(self,model,Wmers_by_seq):
        sum = 0
        #param_gamma = self.param['gamma']
        gamma = model.gamma
        for Wmer_set in Wmers_by_seq:
            print len(Wmer_set)
            param_lambda = float(gamma) / len(Wmer_set[0].src)
            Qi = 0 

            for wmer in Wmer_set:
                sum = sum + wmer.Z * wmer.score * wmer.count
                Qi  = Qi + wmer.Z
                #if wmer.Z > 0.5:
                #    print "Z:%f s:%f c:%f\t%s\t%f"%(wmer.Z,wmer.score,wmer.count,wmer,sum)
            sum = sum + (1 - Qi) * Wmer_set[0].srcQ
            sum = sum + (1 - Qi) * math.log(gamma)/math.log(2)
            sum = sum +      Qi  * math.log(param_lambda)         /math.log(2)

        return(sum)
                

    def EM_C(self, Model, probes, store_Zs=''):
        width = Model.width
        bg    = self.bgprob

        '''Build Probelist'''
        c_Probelist = MDsupport.Probelist()
        for probe in probes:
            c_Probelist.append(probe.c_intA, len(probe), probe.logP, probe.c_wmerbgs)

        '''Build C++ PSSM'''
        c_Model          = MDsupport.SeqMat(width)
        c_Model.gamma    = self.param['gamma']
        c_Model.gammawt  = 0.8
        c_Model.deltamin = self.deltamin
        c_Model.beta     = self.beta
        c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'])
        
        '''Gap maddness'''
        for i in range(width):
            c_Model.setmask(i,self.mask[i])

        '''Copy python PSSM data into C++ PSSM'''
        Ljs = zip(['A','C','G','T'],[0,1,2,3])
        for i in range(width):
            for L,j in Ljs:
                c_Model.set(i,j, Model.logP[i][L])

        '''Do EM until convergence (c_Model.deltamin)'''
        c_Model.EMstep(c_Probelist, self.param['gamma'])

        '''Store the Zs with each Wmer, which is useful for doing'''
        '''some post-processing later.'''
        if 0 and store_Zs:
            print '#Re-organizing Z data...',;sys.stdout.flush()
            for Wmer_set,i in zip(Wmers_by_seq, range(len(Wmers_by_seq))):
                for wmer,j in zip(Wmer_set, range(len(Wmer_set))):
                    wmer.Z = c_Probelist.get_Z(i,j)
            print 'Done.'; sys.stdout.flush()

        '''Build new python Motif from data in updated C++ Motif'''
        newW = []
        Ljs = zip(['A','C','G','T'],[0,1,2,3])
        _t = 0
        for i in range(width):
            d = {}
            for L,j in Ljs:
                d[L] = c_Model.get_c(i,j)
                _t   = _t + d[L]
            newW.append(d)


        if _t == 0:
            return None

        Mtmp = MotifTools.Motif(None,self.bgprob)
        Mtmp.compute_from_counts(newW,0.001)
        Mtmp.gamma = c_Model.gamma
        Mtmp.joint = c_Model.joint

        '''Pass it back'''
        return(Mtmp)

    def best_by_Z(self,Wmers_by_seq):
        ans = []
        for Wmer_set in Wmers_by_seq:
            bestwmer,bestZ = '', 0
            for wmer in Wmer_set:
                if wmer.Z > bestZ:
                    bestZ = wmer.Z
                    bestwmer = wmer.seq
            if bestwmer: ans.append(bestwmer)
            if 0 and self.verbose: print "Best: %s %f"%(bestwmer,math.log(bestZ))
        return(ans)

    def all_Wmers(self,N,seq):
        forw = []
        rev  = []
        seqrc = MotifTools.revcomplement(seq)
        Mlh = theMarkovBackground.highestorder
        Mlb = theMarkovBackground.logbackground
        MCP = theMarkovBackground.CP
        Fbg = Mlb(seq)
        Rbg = Mlb(seqrc)
        nmask = map(lambda x:1-x, self.mask)

        '''
        ?? QUESTION: Is it sensible to compute the background probabilities
        this way?
        
        1) BG of complementary strand is taken as equal to primary strand.
        2) Letters inside the motif window are not used for conditional probabilities.
           As a result, the calculation essentially breaks down to the log probability the
           background emits the sequence to the left of the window plus the log probability
           the background emits the sequence to the right.
        3) I\'ve worked out an efficient way to compute this by
           a) Compute the background probability for the entire probe/sequence
           b) (Quick) Compute logQdiff below
           c) Subtract
        '''

        for i in range(len(seq)-N+1):
            subseq = seq[i:i+N]

            '''Build Wmer information'''
            #Wtmp        = Wmer(subseq)
            left        = seq[0:i]
            right       = seq[i+N:]
            #Wtmp.lflank = left
            #Wtmp.rflank = right
            #if i==0: Wtmp.src    = seq
            #Wtmp.srcQ   = Fbg
            #Wtmp.i      = i

            '''This is the fast way'''
            logQdiff = Mlb(left[-Mlh:] + subseq + right[0:Mlh]) - Mlb(left[-Mlh:]) - Mlb(right[0:Mlh])
            logQtot = Fbg - logQdiff

            '''Add a bit back for intervening bases in the "gap" '''
            gapbg = 0
            for p in range(N):
                gapbg = gapbg + MCP[subseq[p]] * nmask[p]
            logQtot = logQtot + gapbg

            '''Build Wmer-reverse complement information'''
            #Wtmprc = Wmer(Wtmp.rc)
            #Wtmprc.lflank = seqrc[0:-(i+N)]  #Check this in case it is ever necessary
            #if i!=0:
            #    Wtmprc.rflank = seqrc[-i:]   #Necessary [11-12-02]
            #else:
            #    Wtmprc.rflank = ''
            #Wtmprc.logQtot = Wtmp.logQtot
            #Wtmprc.srcQ    = Wtmp.srcQ
            #Wtmprc.i       = i
            forw.append(logQtot)
            rev.append(logQtot)
        W = []
        W.extend(forw)
        W.extend(rev)
        #seq.c_wmerbgs = MDsupport.list2double(map(lambda x: x.logQtot, W))
        #MDsupport.printdouble(seq.c_wmerbgs,len(W))
        return(W)

        
class MarkovBackground:
    def __init__(self,species='YEAST',seqs=''):
        if  seqs:
            self.sourcefile = 'Runtime (%d sequences)'%len(seqs)
        elif  species.find('.6MBG') >= 0:
            self.sourcefile = species
        elif species[0:5] == 'YEAST':
            self.sourcefile = TAMO.paths.Whiteheaddir+'Yeast6kArray/yeast.intergenic.6.freq'
            TAMO.paths.CHECK(self.sourcefile,'Whitehead')
        elif species[0:5] == 'HUMAN':
            self.sourcefile = TAMO.paths.Whiteheaddir+'Human13kArray/human_elongated_probbesQC250.6MBG'
            TAMO.paths.CHECK(self.sourcefile,'Whitehead')
        elif os.path.exists(re.sub('.fsa|.fasta','.6MBG',species)):
            self.sourcefile = re.sub('.fsa|.fasta','.6MBG',species)
        elif os.path.exists(species) and (species.find('.fsa') >=0):
            self.sourcefile = species
            print "EM.MarkovBackground: Computing background from %s"%species
            sys.stdout.flush()
            self.freqs_from_seqs(Fasta.seqs(species))
        #elif os.path.exists(species):
        #    self.sourcefile = species
        else:
            print 'EM.MarkovBackground: Unknown species %s, using Yeast'
            self.sourcefile = TAMO.paths.Whiteheaddir+'Yeast6kArray/yeast.intergenic.6.freq'
            TAMO.paths.CHECK(self.sourcefile,'Whitehead')
        self.species = species
        self.D  = {}
        self.F  = {}         #Frequencies
        self.CP = {}         #log2(Conditional Probabilities)  CP['ACTG'] = p( G | ACT ) 
        self.nmers_by_size = map(lambda x:[],range(0,10))
        self.highestorder = 0
        if seqs:
            print "EM.MarkovBackground: Computing background from %d sequences"%len(seqs)
            self.freq_from_seqs(seqs)
        else:
            self.freq_from_file()
        self.compute_conditional()
        self.totD = {}
        #self.totD = shelve.open('MarkovBackgroundTotals'+`self.highestorder`+'.shelve')

    def zeroth(self):
        D = {}
        for L in ['A', 'C', 'G', 'T']:
            D[L] = self.F[L]
        return(D)

    def permute(self,letters, depth, seqs=[''],curdepth=0):
        newseqs = []
        for seq in seqs:
            for letter in letters:
                newseqs.append(seq + letter)
        if depth > curdepth:
            return(self.permute(letters,depth,newseqs,curdepth + 1))
        else:
            return(seqs)

    def compute_conditional(self):
        for lett in ['A','C','T','G']:
            self.CP[lett]= math.log(self.F[lett]) / math.log(2)
        for depth in range(2,self.highestorder+1):
            for leading_seq in self.nmers_by_size[depth-1]:
                sub_total = self.F[leading_seq]
                sub_tot   = 0
                for trailing_lett in ['A','C','T','G']:
                    sub_key = "%s%s"%(leading_seq,trailing_lett)
                    if self.F.has_key(sub_key):
                        sub_tot = sub_tot + self.F[sub_key]
                for trailing_lett in ['A','C','T','G']:
                    sub_key = "%s%s"%(leading_seq,trailing_lett)
                    if self.F.has_key(sub_key):
                        self.CP[sub_key] = math.log(self.F[sub_key]/sub_tot) / math.log(2)
        
    def freq_from_fasta(self,fastafile):
       seqsD = Fasta.load(sys.argv[1])
       seqs  = seqsD.values()
       self.freq_from_seqs(seqs)

    def freq_from_seqs(self,seqs):
       self.highestorder = 6
       for w in range(1,7):
           allnmers = permute(w)
           nmersT = MotifTools.top_nmers(w,seqs,'with counts','purge Ns')
           self.nmers_by_size[w] = allnmers[:]
           nmersD = {}
           total = 0.0
           for nmer in allnmers: #Pseudo count
               nmersD[nmer] = 1 
               total = total + 1
           for nmer,count in nmersT:
               try: 
                   rc = MotifTools.revcomplement(nmer)
                   nmersD[nmer] = nmersD[nmer] + count
                   nmersD[rc]   = nmersD[rc]   + count
                   total = total + 2*count
               except KeyError:
                   pass
           for nmer in nmersD.keys():
               rc = MotifTools.revcomplement(nmer)
               f  = nmersD[nmer]/total
               self.F[nmer] = f
               self.F[rc]   = f

    def freq_from_seqs_old(self,seqs):
        self.highestorder = 4
        for depth in range(1,6):
            nmersT = MotifTools.top_nmers(depth, seqs, "TUPLES")
            self.nmers_by_size[depth] = map(lambda x:x[0],nmersT)
            total = 0
            for nmer,count in nmersT:
                total = total + count
            for nmer,count in nmersT:
                rc = MotifTools.revcomplement(nmer)
                if nmer == rc:                       #correct top_nmers 
                    f   = float(count)/total         #palindrome count
                else:
                    f   = float(count)/total/2
                self.F[nmer] = f
                self.F[rc]   = f
        for depth in range(0):                       #For debugging
            total = 0
            for k in self.F.keys():
                if len(k) == depth:
                    total = total + self.F[k]
                    print k, self.F[k]
            print depth,total
                
    def study_seqs(self,seqs):
        for depth in range(1,6):
            nmersT = MotifTools.top_nmers(depth, seqs, "TUPLES")
            total = 0
            for nmer,count in nmersT:
                total = total + count
                rc = MotifTools.revcomplement(nmer)
            for nmer,count in nmersT:
                f   = math.log(float(count)/total)/math.log(2)
                f_2 = math.log(0.5 * float(count)/total)/math.log(2)
                rc = MotifTools.revcomplement(nmer)
                if rc != nmer:
                    self.D[nmer] = f_2
                    self.D[rc]   = f_2
                else:
                    self.D[nmer] = f
        for depth in range(0):
            total = 0
            for k in self.D.keys():
                if len(k) == depth:
                    total = total + pow(2,self.D[k])
                    print k, pow(2,self.D[k])
            print depth,total
        self.highestorder = 5

    def freq_from_file(self):
        print "# Loading High-order Markov background model from:\n#\t\t%s"%self.sourcefile
        sys.stdout.flush()
        FH = open(self.sourcefile,'r')
        lines = FH.readlines()
        FH.close()
        for line in lines:
            if line[0] == '#':  #Comment line
                continue
            (key,freq) = line.split()
            self.F[key.upper()] = float(freq)
        self.highestorder = -1 + max(map(len,self.F.keys()))
        for depth in range(1,self.highestorder):
            a = []
            for key in self.F.keys():
                if len(key) == depth:
                    a.append(key)
            #print depth,len(a),self.nmers_by_size
            self.nmers_by_size[depth] = a
    def logbackground(self,seq):
        #temp_order = self.highestorder
        #self.highestorder=6
        _T = 0
        if  len(seq) == 0:
            _T = 0
        elif self.totD.has_key(seq):   #Memo-ize frequent lookups
            _T = self.totD[seq]
        else:
            for i in range(len(seq)):
                start  = max(0,i+1-self.highestorder)
                subseq = seq[start:i+1]
                #print start,i+1-self.highestorder,subseq,self.D[subseq]
                #print "subSEQ: ",start,i+1,subseq
                try:
                    _T = _T + self.CP[subseq]
                except KeyError:
                    print "Error analysing:\n%s"%seq
                    print "<%s>  %d"%(subseq, len(self.CP))
                    print "start=%d   i=%d  ho=%d"%(start,i,self.highestorder)
                    _T = _T + self.CP[subseq]
            self.totD[seq] = _T
        #self.highestorder=temp_order
        return(_T)
    def _keys(self):
        return self.D.keys()
    def ___getitem__(self,key):
        return self.D[key]
    def _has_key(self):
        return D.has_key(key)

class Zeroth(MarkovBackground):
    def __init__(self,bgD={'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25, }):
        self.D  = {}
        self.F  = {}         #Frequencies
        self.CP = {}         #log2(Conditional Probabilities)  CP['ACTG'] = p( G | ACT ) 
        self.nmers_by_size = map(lambda x:[],range(0,10))
        self.highestorder = 1
        self.F = bgD
        self.compute_conditional()
        self.totD = {}

class Probe:
    '''
    Probe object: Extension of string class containing
    extra information, including logP and other information
    '''
    def __init__(self,seq=''):
        self.seq    = seq
        self.info   = ''
        self.logP   = 0
        self.c_wmerbgs= None
        if seq:
            self.c_intA = MDsupport.seq2int(seq)
            self.compute_logP()
    def compute_logP(self):
        global theMarkovBackground
        self.logP  = theMarkovBackground.logbackground(self.seq)
    def __cmp__(self, other):
        if type(other) == type(self):
            seq = other.seq
        else:
            seq = other
        return cmp(self.seq,seq)
    def __repr__(self):
        return self.seq
    def __len__(self):
        return len(self.seq)
    def __getitem__(self,n):
        return self.seq[n]
    def __getslice__(self,i,j):
        return self.seq[i:j]
    def __hash__(self):
        return(hash(self.seq))
    def translate(self,table):
        return self.seq.translate(table)
    def __commented_del__(self):
        print "PROBE: del (len %d)"%len(self.seq)
        
def loadMarkovBackground(species='YEAST',seqs=[]):
    global theMarkovBackground
    theMarkovBackground = MarkovBackground(species,seqs)

def permute(depth, letters=['A','C','G','T'], seqs=[''],curdepth=0):
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(depth,letters,newseqs,curdepth + 1))
    else:
        return(seqs)

def log2_sum(logx, logy):
    #ans = logx + math.log(1 + math.exp(logy - logx))
    if (logx > logy):
        (loga,logb) = (logx,logy)
    else:
        (loga,logb) = (logy,logx)
    if  logb < -1e100:
        ans = loga
    else:
        ans = loga + math.log(1 + math.pow(2,logb - loga))/math.log(2)
    return(ans)

log2_sum = MDsupport.log2_sum    
