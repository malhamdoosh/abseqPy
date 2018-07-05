#!env python
import sys, pickle, random, math, re
from TAMO import MotifTools
from TAMO.MD import EM, MDsupport
from TAMO.seq import Fasta

def main():
    ##########################################################################################
    #THEME.py: THEME module for performing cross-validated hypothesis testing on transcription
    #factor binding data.
    #Usage: python THEME.py foreground_fasta_file (file path) background_fasta_file (file path)
    #hypothesis_index (integer)  -fse hypothesis_file (file path) -markov markov_background (file path)
    #-motif_file output_file (file path) -cv fold cross-validation (integer)
    ##########################################################################################

    if (len(sys.argv)<4):
        print "Usage: THEME.py foreground.fsa background.fsa hypotheses.txt"
        sys.exit(1)

    fg_file = sys.argv[1]           #get fasta file with foreground sequences
    bg_file = sys.argv[2]           #get fasta file with background sequences
    test_indices = sys.argv[3]      #colon separated indices into fse file
    cv_level = 2                    #default 2-fold cross-validation
    refine = 1
    randomize = 0
    beta = 0.0
    delta = 0.001
    motif_file = 'dummy.out'
    dump_categories_to_file = 0
    test_family = ''
    
    #read in any command line options
    for arg, i in zip(sys.argv,range(len(sys.argv))):
        if (arg == '-cv'):
            cv_level = int(sys.argv[i+1])
        if (arg == '-markov'):
            markov_file = sys.argv[i+1]
        if (arg == '-fse'):
            fse_file = sys.argv[i+1]
        if (arg == '-norefine'):
            refine = 0
        if (arg == '-beta'):
            beta = float(sys.argv[i+1])
        if (arg == '-delta'):
            delta = float(sys.argv[i+1])
        if (arg == '-randomization'):
            randomize = 1
        if (arg == '-motif_file'):
            motif_file = sys.argv[i+1]
        if (arg == '-dump'):
            dump_categories_to_file = 1
        if (arg == '-family'):
            test_family = family
    FH = open(motif_file, 'w')
    FH.write("******THEME Motif Output******")
    FH.close()
    
    random.seed()

    cross_val = THEME(fg_file, bg_file, cv_level, markov_file)
    if ((beta>0.0)and(beta<1.0)) : cross_val.beta = beta/(1-beta)
    cross_val.delta = delta
    cross_val.refine = refine
    cross_val.randomize = randomize
    cross_val.motif_file = motif_file
    if (test_family): cross_val.family = test_family
    if (dump_categories_to_file):
        cross_val.dump = 1

    ###################################################################################
    #get seed sequences that will be tested
    ###################################################################################
    models = []
    fses = MotifTools.load(fse_file)
    if (test_indices=='all'):
        indices = range(len(fses))
    else:
        indices = []
        ivals = test_indices.split(':')
        for v in ivals:
            indices.append(int(v))
    for i in indices:
        ll = fses[i].logP
        bg = EM.theMarkovBackground.zeroth()
        for pos in ll:
            for letter in pos.keys():
                pos[letter] = pos[letter] - math.log(bg[letter])/math.log(2.0)
        adj_bg_model = MotifTools.Motif_from_ll(ll)
        adj_bg_model.source = fses[i].source
        models.append(adj_bg_model)
        
    (m, err) = cross_val.run_CV(models)
        
class THEME:
    def __init__(self, fg_file, bg_file, cv_level, markov_file):
        self.cv_level = cv_level
        self.randomize = 0
        self.beta = 0.0
        self.delta = 0.001
        self.refine = 1
        self.motif_file = 'dummy.out'
        self.dump = 0
        self.family = ''
        self.datafiles = (fg_file,bg_file)
        
        MAX_FG = 2000
        
        #LOAD MARKOV BACKGROUND#
        print "Loading Markov background file from %s"%markov_file
        EM.loadMarkovBackground(markov_file)    

        ##################################################################################
        #divide input sequences into groups according to the desired cross-validation level
        ###################################################################################
        print "Processing input sequences...."
        self.fg_seqs = Fasta.load(fg_file)   #load foreground sequences
        for key in self.fg_seqs.keys():
            fseq = self.fg_seqs[key]
            self.fg_seqs[key] = fseq.split()[0]
        self.all_probes = Fasta.load(bg_file)   #load background sequences
        Fasta.delN(self.fg_seqs)
        Fasta.delN(self.all_probes)

        #first delete any sequences from background that are present in foreground
        for key in self.fg_seqs.keys():
            if (self.all_probes.has_key(key)):
                del self.all_probes[key]

        for key in self.all_probes.keys():
            if ((len(self.all_probes[key])==0) or (re.search('[SWMKRY]', self.all_probes[key]))):
                del self.all_probes[key]
                print "deleting %s"%key
                
        while (len(self.fg_seqs.keys())>MAX_FG):
            del self.fg_seqs[self.fg_seqs.keys()[random.randint(0,(len(self.fg_seqs.keys())-1))]]

    def run_CV(self,models):
        num_bg = len(self.all_probes.keys())
        num_fg = len(self.fg_seqs.keys())
       
        for key in self.fg_seqs.keys():
            self.all_probes[key] = self.fg_seqs[key]

        self.probes = self.all_probes.keys()
                
        if (self.randomize): trials = 50
        else: trials = 1

        self.models = models
        
        std = 0.0
        mean = 0.0
        sum_sq = 0.0
        sum_mean = 0.0
        for trial in range(trials):
            if (self.randomize):
                N = num_bg/num_fg
                if (N>10): N = 10
                fg = []
                bg = []
                total = len(self.probes)
                while(len(fg)<num_fg):
                    sp = random.randint(0,(total-1))
                    if (not(fg.count(self.probes[sp])>0)):
                        fg.append(self.probes[sp])
                while(len(bg)<(N*num_fg)):
                    sp = random.randint(0,(total-1))
                    if ( (not(bg.count(self.probes[sp])>0)) and (not(fg.count(self.probes[sp])>0)) ):
                        bg.append(self.probes[sp])
            else:            
                #select a random under-sampled set of background sequences            
                N = num_bg/num_fg
                (N, bg) = self.undersample(N)
                num_bg = N*num_fg
                fg = self.fg_seqs.keys()

            fg_group_size = len(fg)/self.cv_level   #determine foreground group size
            bg_group_size = len(bg)/self.cv_level   #determine background group size

            #separate sequences into groups
            fg_groups = {}
            bg_groups = {}
            temp_fg = []
            temp_bg = []
            for a in fg: temp_fg.append(a)
            for b in bg: temp_bg.append(b)
            for i in range(self.cv_level):
                fg_groups[i] = []
                bg_groups[i] = []
                if (i==(self.cv_level-1)):
                    fg_groups[i] = temp_fg[0:]
                    bg_groups[i] = temp_bg[0:]
                else:
                    for j in range(fg_group_size):
                        entry = temp_fg[random.randint(0,len(temp_fg)-1)]
                        temp_fg.remove(entry)
                        fg_groups[i].append(entry)
                    for j in range(bg_group_size):
                        entry = temp_bg[random.randint(0,len(temp_bg)-1)]
                        temp_bg.remove(entry)
                        bg_groups[i].append(entry)

        ####################################################################################
        #for each seed we run EM to train a probability model, calculate the maximum LLR for
        #each input sequence in a group, train a SVM classifier on the training data, and
        #calculate classification error on the test set.  We repeat this for each cv-group
        #and determine the mean cross-validation error for each hypothesis.
        ####################################################################################
            self.classification_errors = {}
            self.refined_motifs = {}
            classifier = {}
            self.best_motif = []
            num_models = len(self.models)
            for k in range(num_models):
                self.best_motif.append(None)
                (self.refined_motifs[k], self.classification_errors[k]) = self.train_classifier(k, fg_groups, bg_groups, N, trial)
               
                #print out some information
                best_c = 0
                best_beta = 0.0
                min = 1.0
                mean = 0.0
                for beta in self.classification_errors[k].keys():
                    for j in range(len(self.classification_errors[k][beta][0])):
                        mean = 0.0
                        for i in range(self.cv_level):
                            mean = mean + self.classification_errors[k][beta][i][j]
                        mean = mean/self.cv_level
                        if (mean<min):
                            min = mean
                            best_c = j
                            best_beta = beta
                if (self.randomize):
                    sum_sq = sum_sq + (min*min)
                    sum_mean = sum_mean + min
                    cv_mean = sum_mean/(trial+1)
                    if (trial>0):
                        try: std = math.sqrt((sum_sq - (sum_mean*sum_mean)/(trial+1))/(trial))
                        except: std = 0.0
                    stddev_err = 0.71*std/math.sqrt(trial+1)
                else:
                    (self.best_motif[k], classifier[k], fn) = self.train_final(k, fg, bg, N, best_beta)
                    print "\r\r---------------New Hypothesis---------------"
                    print "Unrefined Hypothesis: %s"%(self.models[k])
                    print "Refined Hypothesis: %s Optimal Beta: %f"%(self.best_motif[k], best_beta/(best_beta+1.0))
                    print "Mean %i-fold cv error: %f"%(self.cv_level,min)
                    LLR_thresh = (classifier[k][2]*classifier[k][1]/classifier[k][0] - classifier[k][3])
                    print "LLR match threshold: %f True positives: %i False Negatives: %i"%(LLR_thresh/self.best_motif[k].maxscore, (len(fg)-fn), fn)
                    self.best_motif[k].source = self.models[k].source
                    if (self.family!=''): self.best_motif[k].family = self.family
                    self.best_motif[k].dataset = self.datafiles[0]
                    self.best_motif[k].bgfile = self.datafiles[1]
                    self.best_motif[k].beta = best_beta/(best_beta+1.0)
                    self.best_motif[k].match_thresh = LLR_thresh
                    self.best_motif[k].cverror = min
                    print "Motif matches in positive input set:"
                    best_pssm = MDsupport.Motif2c_PSSM(self.best_motif[k])
                    for seq in fg:
                        sites = []
                        matches = best_pssm.matchstarts(self.all_probes[seq],LLR_thresh)
                        if (matches):
                            line = seq + '------>  '
                            for match in matches:
                                entry = str(match) + ': ' + self.all_probes[seq][match:match+self.best_motif[k].width] + ' '
                                line = line + entry
                            print line
                        
        if (self.randomize):
            print "Mean: %f Std.Dev: %f Error: %f Percent error: %f"%(cv_mean,std,stddev_err,(stddev_err/std))
        else:
            MotifTools.save_motifs(self.best_motif, self.motif_file)
        return((self.best_motif, min))

    def undersample(self,N):
        num_fg = len(self.fg_seqs.keys())
        seqs = self.all_probes.keys()
        if (N>10):
            N = 10
            temp_bg = []
            for i in range(N*num_fg):
                b_key = self.fg_seqs.keys()[0]
                while (self.fg_seqs.has_key(b_key)):
                    s_index = random.randint(0, len(seqs)-1)
                    b_key = seqs[s_index]
                temp_bg.append(b_key)
                seqs.remove(b_key)
        else:
            num_bg = len(self.all_probes.keys())
            for i in range((num_bg-N*num_fg)):
                b_key = self.fg_seqs.keys()[0]
                while (self.fg_seqs.has_key(b_key)):
                    s_index = random.randint(0, len(seqs)-1)
                    b_key = seqs[s_index]
                seqs.remove(b_key)
            temp_bg = seqs
        return((N,temp_bg))

    def train_classifier(self, model_num, fg_groups, bg_groups, N, trial):
        refined_models = {}
        errors = {}
        if (self.refine):
            if (self.beta > 0.0): betas = [self.beta]
            else: betas = [(0.05/0.95), (0.1/0.9), 0.5, 1.0, 2.0, 1000000000000000000]
        else: betas = [1000000000000000000]
        for beta in betas:
            refined_models[beta] = []
            errors[beta] = []
            for i in range(self.cv_level):
                train_fg = []
                train_bg = []

                #separate into training and test sets
                for j in range(self.cv_level):
                    if (j!=i):
                        for s in fg_groups[j]:
                            train_fg.append(s)
                        for s in bg_groups[j]:
                            train_bg.append(s)
                test_fg = fg_groups[i]   #positive test set
                test_bg = bg_groups[i]   #negative test set

                #get input sequences for this run
                input_seqs = []
                for s in train_fg:
                    iseq = self.all_probes[s].upper()
                    iseq = re.sub(";","",iseq)
                    if (re.search("N", iseq)):
                        iseq = re.sub("N","",iseq)
                    if (len(iseq)>0): input_seqs.append(iseq)

                #print input_seqs
                #Run TAMO-EM using input sequences and current seed
                if (self.refine):
                    refined_models[beta].append(self.train_model(model_num, input_seqs, beta))
                else:
                    refined_models[beta].append(self.models[model_num])
                #Calculate maximum log-likelihood ratio for each sequence
                train_pos = self.get_LLRs(refined_models[beta][i], train_fg)
                train_neg = self.get_LLRs(refined_models[beta][i], train_bg)
                test_pos  = self.get_LLRs(refined_models[beta][i], test_fg)
                test_neg  = self.get_LLRs(refined_models[beta][i], test_bg)

                #perform SMOTE over-sampling of positive data set
                over_sampled_positive = self.SMOTE([train_pos, test_pos], N, N)
                train_pos = over_sampled_positive[0]
                test_pos = over_sampled_positive[1]
                
                the_mean = 0.0
                for val in train_pos:
                    the_mean = the_mean + val
                #print "pos: %f"%(the_mean/len(train_pos))
                the_mean = 0.0
                for val in train_neg:
                    the_mean = the_mean + val
                #print "neg: %f"%(the_mean/len(train_neg))

                #Train SVM to classify our training set
                c_vals = [1.0e-10, 1.0e-4, 1.0e-3, 1.0e-2, 0.05, 0.1, 1.0, 10.0, 100.0]
                cv_errs = []
                for c in c_vals:
                    classifier = self.SVM_train(train_pos, train_neg, c)
                    train_err = self.SVM_test(classifier, train_pos, train_neg)
                    #print "cv: %i c-val: %f training error: %f"%(i, c, train_err)

                    # Use the newly trained classifier to classify the test set
                    cv_errs.append(self.SVM_test(classifier, test_pos, test_neg))
                errors[beta].append(cv_errs)
        return((refined_models, errors))

    def train_final(self, model, fg, bg, N, beta):
        input_seqs = []
        for s in fg:
            iseq = self.all_probes[s].upper()
            iseq = re.sub(";","",iseq)
            if (re.search("N", iseq)):
                iseq = re.sub("N","",iseq)
            if (len(iseq)>0): input_seqs.append(iseq)

        if (self.refine):
            final_motif = self.train_model(model, input_seqs, beta)        
        else:
            final_motif = self.models[model]
        train_pos = self.get_LLRs(final_motif, fg)
        train_neg = self.get_LLRs(final_motif, bg)
        over_sampled_positive = self.SMOTE([train_pos], N, N)[0]

        #Train SVM to classify our training set
        c_vals = [1.0e-10, 1.0e-4, 1.0e-3, 1.0e-2, 0.05, 0.1, 1.0, 10.0, 100.0]
        best_classifier = None        
        best_err = 1.0
        for c in c_vals:
            classifier = self.SVM_train(over_sampled_positive, train_neg, c)
            train_err = self.SVM_test(classifier, over_sampled_positive, train_neg)
            if (train_err<best_err):
                best_err = train_err
                best_classifier = classifier
        (train_err, fp, fn) = self.SVM_test(best_classifier, train_pos, [], 1)
        if (self.dump):
            motif = {}
            no_motif = {}
            for name, val in zip(fg,train_pos):
                train_err = self.SVM_test(best_classifier, [val], [])
                if (train_err):
                    no_motif[name] = self.all_probes[name]
                else:
                    motif[name] = self.all_probes[name]
            motif_fsa = self.motif_file.split('.')[0] + '.pos.fsa'
            no_motif_fsa = self.motif_file.split('.')[0] + '.neg.fsa'
            Fasta.write(motif, motif_fsa)
            Fasta.write(no_motif, no_motif_fsa)
        return((final_motif, best_classifier, fn))

    def train_model(self, hypothesis_num, inputs, beta=0.01):
        if (beta==1000000000000000000):
            return(self.models[hypothesis_num])
        else:
            theEM = EM.EM([],[],None,"")
            theEM.models = [self.models[hypothesis_num]]
            theEM.beta = beta 
            theEM.deltamin = self.delta
            theEM.param['gamma'] = 0.75
            theEM.seqs.extend(inputs)
            theEM.EM_Cstart()
            return(theEM.candidates[0].pssm)

    #################################################################################################
    # method to evaluate the maximum LLR score for a set of sequences
    #################################################################################################
    def get_LLRs(self, motif, seq_tags, verbose=0):
        LLRs = []
        pssm = MDsupport.Motif2c_PSSM(motif)
        for key in seq_tags:
            LLRs.append(pssm.scanbest(self.all_probes[key]))
            '''
            seq_fwd = self.all_probes[key]                         #get the sequence data
            seq_rev = str(MotifTools.revcomplement(seq_fwd))
            try:
                score_fwd = pssm.scanbest(seq_fwd)
                score_rev = pssm.scanbest(seq_rev)
            except:
                print seq_fwd
                print seq_rev
            LLRs.append(max(score_fwd,score_rev))
            '''
        return(LLRs)

    ###################################################################################################
    #method to upsample data using the SMOTE algorithm
    ###################################################################################################
    def SMOTE(self, data, N, k):

        #assemble data into a sorted list labeled by group
        data_list = []
        for i in range(len(data)):
            for point in data[i]:
                data_list.append((point, i))
        data_list.sort()

        #generate list of k nearest neighbours for each point
        nearest_neighbours = {}
        for i in range(len(data_list)):
            nearest_neighbours[i] = []
            distances = []
            for j in range(k): distances.append((10000.0, 10000.0))
        
            start = i-k
            if (start<0): start = 0
            end = i+k+1
            if (end>=len(data_list)): end = len(data_list)

            for j in range(start,end):
                if (i!=j):
                    distance = data_list[j][0]-data_list[i][0]
                    abs_dist = math.fabs(distance)
                    rem = 0
                    for d in distances:
                        if (abs_dist<d[0]):
                            rem = d
                            break
                    if (rem!=0):
                        distances.remove(rem)
                        distances.append((abs_dist, distance))
                        distances.sort()
                        distances.reverse()

            nearest_neighbours[i] = distances

        #generate upsampled data
        up_data = []
        for i in range(len(data)):
            up_data.append([])
        for i in range(len(data_list)):
            up_data[data_list[i][1]].append(data_list[i][0])
            for j in range(N-1):
                #generate random integer
                neighbour = random.randint(0,(k-1))
            
                #generate new data point
                new_val = data_list[i][0] + nearest_neighbours[i][neighbour][1]*random.random()

                up_data[data_list[i][1]].append(new_val)
        return(up_data)

    ###################################################################################################
    #method to train a 1-D SVM classifier
    ###################################################################################################
    def SVM_train(self, positives, negatives, _c):

        #scale the data
        c = _c
        maxval = max(max(positives),max(negatives))
        minval = min(min(positives),min(negatives))
        den = (maxval-minval)/2.0
        bias = (minval-maxval)/2.0 - minval
        ptrain = []
        ntrain = []
        for val in positives:
            ptrain.append(-(val+bias)/den)
        for val in negatives:
            ntrain.append(-(val+bias)/den)

        #sort positives and negatives
        ptrain.sort()
        ptrain.reverse()
        ntrain.sort()
        #ntrain.reverse()

        #calculate cumulative sums
        d_pos = [0.0]
        d_neg = [0.0]
        d = 0.0
        for p in ptrain:
            d = d + p
            d_pos.append(d)
        d_pos.reverse()
        d_pos = d_pos[0:(len(d_pos)-1)]
        d = 0.0
        for n in ntrain:
            d = d + n
            d_neg.append(d)
        d_neg.reverse()
        d_neg = d_neg[0:(len(d_neg))-1]

        ptrain.reverse()
        ntrain.reverse()
        min_L = 1000000000
        for i in range(min(len(ptrain),len(ntrain))):
            p_index = len(ptrain) - 1 - i
            n_index = len(ntrain) - 1 - i
            if (ptrain[p_index]==ntrain[n_index]):
                w = 1.0e6
            else:
                w = 2.0/(ptrain[p_index]-ntrain[n_index])
            b = (1.0-ptrain[p_index])/(w*w)
            sum_etas = 2*i - w*(d_pos[p_index] - d_neg[n_index])
            L = (w*w)/2.0 + c*sum_etas
            if (L<min_L):
                min_L = L
                best_w = w
                best_b = b
        best_classifier = (best_w, best_b, den, bias)
        return(best_classifier)

    ####################################################################################################
    #method to test a classifier on some data
    ####################################################################################################
    def SVM_test(self, classifier, positives, negatives, detailed = 0):
        w = classifier[0]       #slope of our classifier
        b = classifier[1]       #y-intercept of our classifier
        scale = classifier[2]   #for scaling inputs
        offset = classifier[3]  #for scaling inputs
        fp = 0                  #false positives
        fn = 0                  #false negatives

        for val in positives:
            label = w*(val+offset)/scale - b
            if (label<0.0): fn = fn + 1
        for val in negatives:
            label = w*(val+offset)/scale - b
            if (label>0.0): fp = fp + 1

        err = float((fp+fn))/(len(positives)+len(negatives))
        #print "false positives: %i false negatives: %i error: %f"%(fp,fn,err)
        if (detailed):
            return(err, fp, fn)
        else:
            return(err)

if __name__ == '__main__': main()
