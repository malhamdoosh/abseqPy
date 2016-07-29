from igRepUtils import writeCountsToFile, \
    fastq2fasta, mergeReads, writeListToFile
from collections import Counter
from os.path import exists
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from pandas.core.frame import DataFrame
from pandas.io.parsers import read_csv
from NGSmotifs import generateMotifs
from numpy import Inf, random, isnan, NaN, logical_not
from igRepUtils import  writeCountsCategoriesToFile
import gc
from igRepUtils import compressSeqGeneLevel, calMaxIUPACAlignScores, compressSeqFamilyLevel, loadIGVSeqsFromFasta, findBestMatchedPattern
from igRepUtils import compressCountsFamilyLevel, findHitsRegion, replaceIUPACLetters, extractProteinFrag, findBestAlignment,  findHits
from NGSmotifs import generateProteinLogos
from Bio.SeqRecord import SeqRecord
from pandas.core.series import Series
from pandas.tools.merge import concat
from igRepUtils import compressCountsGeneLevel
from config import MEM_GB
import traceback
from igRepAuxiliary import annotateIGSeqRead, writeAbundanceToFiles
from igRepPlots import plotSeqLenDist, plotSeqLenDistClasses, plotVenn, plotDist,\
    plotSeqDuplication, plotSeqDiversity


FR4_CONSENSUS = "WGQGTLVTVSS"

class IgRepertoire:    
    def __init__(self, args):        
        self.outputDir = args['o']
        self.threads = args['threads']
        self.primer = args['primer']
        self.db = args['db']
        self.bitScore = args['bitscore']
        self.alignLen = args['alignlen']
        self.sStart = args['sstart']
        self.seqType = args['seqtype']
        self.format = args['fmt']
        self.chain = args['chain']
        self.name = args['name']
        if (args['task'] in ['secretion', '5utr']):
            self.upstream = args['upstream']
        if (args['task'] in ['enzymes', 'enzymesimple']):
            self.sitesFile = args['sites']
        if (args['task'] == 'diversity'):
            self.end5 = args['5end']
            self.end3 = args['3end']
            self.actualQstart = args['actualqstart']            
            self.end5offset = args['5endoffset']
            self.fr4cut = args['fr4cut']
        self.readFile1 = args['f1']
        self.readFile2 = args['f2']                                
        self.merge = args.get('merge')
        self.merger = args.get('merger')
        
        self.seqsPerFile = 10.0 ** 5  / 2
        self.cloneAnnot = None
        
    def runFastqc(self):
        pass
        
    def identifyClones(self, outDirFilter= None):    
        outDir = self.outputDir + "annot/"
        if (not os.path.isdir(outDir)):
            os.system("mkdir " + outDir)
        self.cloneAnnotFile = outDir + self.name
        self.cloneAnnotFile += "_clones_annot.tab"  
        
        if (exists(self.cloneAnnotFile) and self.cloneAnnot is None):       
            print("Clones annotation file found and is being loaded ... " + 
                  self.cloneAnnotFile.split('/')[-1])
            sys.stdout.flush()
            self.cloneAnnot = read_csv(self.cloneAnnotFile, sep='\t',
                                       header=0, index_col=0)             
        else:                    
            if self.merge != 'yes':            
                self.readFile = self.readFile1                                          
            else:            
                mergedFastq = mergeReads(self.readFile1 , self.readFile2,
                                         self.threads, self.merger, self.outputDir)
                self.readFile = mergedFastq  
            if (not exists(self.readFile)):
                raise Exception(self.readFile + " does not exist!")            
                     
            # Convert FASTQ file into FASTA format
            if self.format == 'fastq':        
                self.readFasta = fastq2fasta(self.readFile, self.outputDir)                
            elif self.format == 'fasta':
                self.readFasta = self.readFile                
            else:
                raise Exception('unknown file format! ' + self.format)            
            sys.stdout.flush()
            # Estimate the IGV family abundance for each library        
            (self.cloneAnnot, filteredIDs) = annotateIGSeqRead(self, self.readFasta,                                                                                                            
                                                          self.seqType)
            sys.stdout.flush()            
            gc.collect()
            writeListToFile(filteredIDs, outDir + self.name + "_unmapped_clones.txt")
            # export the CDR/FR annotation to a file
            print("Clones annotation file is being written to " + self.cloneAnnotFile)         
            self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
        if outDirFilter:    
            ## Filter clones based on bitscore, alignLen and sStart
            print("Clones are filtered based on the following criteria: ")
            print("\tBit score: " + `self.bitScore` )
            print("\tAlignment length: " + `self.alignLen` )
            print("\tV gene start: " + `self.sStart`)
            selectedRows = ((self.cloneAnnot['bitscore'] >= self.bitScore[0]) &  # check bit-Score
                    (self.cloneAnnot['bitscore'] <= self.bitScore[1]) &
                    (self.cloneAnnot['alignlen'] >= self.alignLen[0]) & # check alignment length
                    (self.cloneAnnot['alignlen'] <= self.alignLen[1]) &
                    (self.cloneAnnot['vstart'] >= self.sStart[0]) & # check subject (V gene) start position
                    (self.cloneAnnot['vstart'] <= self.sStart[1]))
            filteredIDs = self.cloneAnnot[logical_not(selectedRows)]            
            filteredIDs = filteredIDs.index.tolist() 
            writeListToFile(filteredIDs, outDirFilter + self.name + "_filtered_out_clones.txt")            
            self.cloneAnnot = self.cloneAnnot[selectedRows]        

    def analyzeAbundance(self):    
        # Estimate the IGV family abundance for each library        
        outDir = self.outputDir + "abundance/"
        if (not os.path.isdir(outDir)):
            os.system("mkdir " + outDir)              
        self.identifyClones(outDir)                 
              
        writeAbundanceToFiles(self.cloneAnnot, self.name, outDir, self.chain)        
        gc.collect()
        
    def analyzeProductivity(self):
        self.identifyClones(True)
        
        self.cdrSeqFile = self.outputDir + "cdrs/" + self.name
        self.cdrSeqFile += "_clones_seq.tab"     
        
        if (not exists(self.cloneAnnotFile)):       
            self.analyzeAbundance('cdrinfo')
            # export the CDR/FR annotation to a file
            self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
            print("Text file has been written to " + self.cloneAnnotFile) 
            # Refine the annotation of CDR3 and FR4
            self.refineCDRInfoAnnotation()  
            self.refineInFramePrediction()
            
        else:
            print("\tCDRs information file was found! ... " + self.cloneAnnotFile.split('/')[-1])
            sys.stdout.flush()
            self.cloneAnnot = read_csv(self.cloneAnnotFile, sep='\t',
                                       header=0, index_col=0)            
            if (not exists(self.cdrSeqFile)):                
                self.refineCDRInfoAnnotation() 
                self.refineInFramePrediction() 
            else:
                print("\tCDRs sequence file was found! ... " + self.cdrSeqFile.split('/')[-1])
                sys.stdout.flush()
                self.cdrSeqs = read_csv(self.cdrSeqFile, sep='\t',
                                       header=0, index_col=0)  
        
    def analyzePrimerSpecificity(self):
        pass
        
    

    def extractUpstreamSeqs(self): 
        if (not exists(self.alignInfoFile)):       
            self.analyzeAbundance('aligninfo')
            self.alignInfo.to_csv(self.alignInfoFile, sep='\t', header=True, index=True)    
            print("Text file has been written to " + self.alignInfoFile)
        else:
            print("\tAlignment information file was found! ... " + self.alignInfoFile.split('/')[-1])
            self.alignInfo = read_csv(self.alignInfoFile, sep='\t',
                                       header=0, index_col=0)
        # extract the upstream DNA sequences and write them into a fasta file  
        print("\tExtracting the upstream sequences ... ") 
        records = []
        revAlign = 0
        trimmedBegin = 0
        expectLength = self.upstream[1] - self.upstream[0] + 1
        trimmedUpstream = 0
        noSeq = 0
        queryIds = self.alignInfo.index
        procSeqs = 0  # processed sequences
        fileHandle = open(self.upstreamFile, 'w')
        fileHandle.close()
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        else:
            records = SeqIO.index(self.readFile1, self.format)
        for id in queryIds:
            record = records[id]  
            
            qsRec = self.alignInfo.loc[record.id]
            if (qsRec.strand != 'forward'):
                revAlign += 1
                record.seq = record.seq.reverse_complement()
            if (qsRec.sstart <= 3):
                end = qsRec.qstart - self.upstream[0] - qsRec.sstart + 1
                if end <= 1:                            
                    noSeq += 1
                else:
                    start = qsRec.qstart - self.upstream[1] - qsRec.sstart + 1
                    if start < 1:
                        start = 1
                    record.description = ""
#                             print(start, end)
                    record.seq = record.seq[int(start - 1):int(end)]
                    if (expectLength != Inf and len(record.seq) < expectLength):
                        trimmedUpstream += 1
                    record.id = record.id + '|' + qsRec.subjid
                    records.append(record)
                    procSeqs += 1
                    if procSeqs % self.seqsPerFile == 0:
                        print('%d/%d sequences have been processed ... ' % (procSeqs, len(queryIds)))
                        SeqIO.write(records, open(self.upstreamFile, 'a'), 'fasta')
                        records = []
            else:
                trimmedBegin += 1
#                         print("The query sequence is not aligned at the start of the IGV sequence! " + record.id)                        
                
                       
        if (len(records) > 0):
            print('%d/%d sequences have been processed ... ' % (procSeqs, len(queryIds)))
            SeqIO.write(records, open(self.upstreamFile, 'a'), 'fasta')     
        if (revAlign > 0):
            print("\t\t\tReversed alignment is not supported ... %d found and excluded!" % (revAlign))
        if (trimmedBegin > 0):
            print("\t\t\tThe query sequence is not aligned within 3bp of the IGV start position ... %d found and excluded!" % (trimmedBegin))
        if (trimmedUpstream > 0):
            print("\t\t\tUpstream sequences shorter than the expected length are detected ... %d found" % (trimmedUpstream))
        if (noSeq > 0):
            print("\t\t\tNo upstream sequence can be extracted, too short, for %d sequences." % (noSeq))
        gc.collect()
            
    def analyzeSecretionSignal(self):        
        print("The diversity of the upstream of IGV genes is being analyzed ... ")
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        self.upstreamFile = self.outputDir + self.name  
        self.upstreamFile += "_upigv_%.0f_%.0f.fasta" % (self.upstream[0],
                                                        self.upstream[1])
        self.alignInfoFile = self.outputDir + self.name
        self.alignInfoFile += "_align_info.tab"
        # extract upstream sequences 
        if (not exists(self.upstreamFile)):
            self.extractUpstreamSeqs()
        else:
            print("\tUpstream sequences file was found! ... " + self.upstreamFile.split('/')[-1])
        self.upstreamFile = os.path.abspath(self.upstreamFile)
        self.alignInfoFile = os.path.abspath(self.alignInfoFile)        
        # plot the distribution of sequence length
        expectLength = self.upstream[1] - self.upstream[0] + 1
        outputFile = self.upstreamFile.replace('.fasta', '_dist.png')
        plotSeqLenDist(self.upstreamFile, self.name, outputFile, self.format)       
        outputFile = self.upstreamFile.replace('.fasta', '_dist_short.png')
        plotSeqLenDist(self.upstreamFile, self.name, outputFile, self.format,
                       expectLength - 1)
        outputFile = self.upstreamFile.replace('.fasta', '_dist_class.png')
        plotSeqLenDistClasses(self.upstreamFile, self.name, outputFile,
                              self.format)
        outputFile = self.upstreamFile.replace('.fasta', '_dist_short_class.png')
        plotSeqLenDistClasses(self.upstreamFile, self.name, outputFile,
                              self.format, expectLength - 1)
        # classify secretion signals based on length, ATG location, gene and gene family
        # analyze intact secretion signals
        self.analyzeSequences(self.name, [expectLength, expectLength], True)
        # analyze trimmed secretion signals 
        self.analyzeSequences(self.name, [1, expectLength - 1], True)
        # analyze 
        
        
    def loadValidSequences(self, sampleName, expectLength, startCodon=True, type='secsig'):
        print("\tSequences between %d and %d are being extracted ... "
              % (expectLength[0], expectLength[1]))
        ighvSignals = {}
        ighvSignalsCounts = Counter()
        ighvSignalsNoATG = []
        noStartCodonCounts = Counter()
        faultyTrans = []
        faultyTransCounts = Counter()
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(self.upstreamFile, 'fasta'))
        else:
            records = SeqIO.index(self.upstreamFile, 'fasta')
        for id in records:
            rec = records[id]
            ighv = rec.id.split('|')[1]
            seq = rec.seq
            if (len(rec) >= expectLength[0] and len(rec) <= expectLength[1]):
                if (not startCodon or "ATG" in seq):
                    if (faultyTransCounts.get(ighv, None) is None):                                                
                        faultyTransCounts[ighv] = 0       
                    if (type == 'secsig'):                 
                        seq = seq.translate(to_stop=False)[1:]
                    if ('X' in seq or '*' in seq):
#                         print(rec.id, str(rec.seq), str(seq))
                        faultyTrans.append(rec)
                        faultyTransCounts[ighv] += 1 
                    elif ('N' not in rec.seq):
#                         rec.seq = seq     
                        if (ighvSignals.get(ighv, None) is None):
                            ighvSignals[ighv] = []    
                            ighvSignalsCounts[ighv] = 0                 
                        ighvSignals[ighv].append(str(rec.seq))  # rec
                        ighvSignalsCounts[ighv] += 1                       
                    else:
                        print('Ignored: ' + str(rec.seq) + ' ' + str(seq))
                        if (type == 'secsig'): 
                            faultyTrans.append(rec)
                            faultyTransCounts[ighv] += 1 
                elif startCodon:                   
                    ighvSignalsNoATG.append(rec)  # seq
                    if (noStartCodonCounts.get(ighv, None) is None):
                        noStartCodonCounts[ighv] = 0
                    noStartCodonCounts[ighv] += 1
        if (sum(ighvSignalsCounts.values()) > 0):
            print("\tThere are %d VALID secretion signals within expected length %s and startCodon=%s " % 
                  (sum(ighvSignalsCounts.values()), str(expectLength), startCodon))
            if (type == 'secsig'):
                title = 'Valid Secretion Signals'
            else:
                title = 'Valid 5`-UTRs'
            writeCountsCategoriesToFile(ighvSignalsCounts, sampleName,
                        self.outputDir + sampleName + '_%s%d%d_valid_' % (type, expectLength[0], expectLength[1]),
                        title)
        # # Faulty secretion signals: stop codons or low quality sequencing
        if (len(faultyTrans) > 0):
            # variant level
            faultySeqFile = self.outputDir + sampleName + '_%s%d%d_faulty_trans.fasta' % (type, expectLength[0], expectLength[1])
            SeqIO.write(faultyTrans, faultySeqFile, 'fasta')            
            writeCountsCategoriesToFile(faultyTransCounts, sampleName,
                        self.outputDir + sampleName + '_%s%d%d_faulty_' % (type, expectLength[0], expectLength[1]),
                        'Faulty Translations')
            print("\tTotal faulty secretion signals is %d (excluded)" % (len(faultyTrans)))
            examples = random.choice(range(len(faultyTrans)), 5)
            for i in examples:
                print(faultyTrans[i].seq, faultyTrans[i].seq.translate())
            faultyTrans = None      
        else:
            faultySeqFile = None    
        # secretion signals with no start codons 
        if (len(ighvSignalsNoATG) > 0):
            noStartCodonFile = self.outputDir + sampleName + '_%s%d%d_no_atg.fasta' % (type, expectLength[0], expectLength[1])
            SeqIO.write(ighvSignalsNoATG, noStartCodonFile, 'fasta')            
            writeCountsCategoriesToFile(noStartCodonCounts, sampleName,
                    self.outputDir + sampleName + '_%s%d%d_no_atg_' % (type, expectLength[0], expectLength[1]),
                    'Secretion Signals without Start Codon')            
            print("\tThere is no ATG codon in %d sequences (excluded). " % (len(ighvSignalsNoATG)))
            examples = random.choice(range(len(ighvSignalsNoATG)), 5)
            for i in examples:
                print(ighvSignalsNoATG[i].seq)
            ighvSignalsNoATG = None
        else:
            noStartCodonFile = None
        gc.collect()
        return (ighvSignals, faultySeqFile, noStartCodonFile)
        
    def analyzeSequences(self, sampleName, expectLength, startCodon=True,
                         type='secsig', clusterMotifs=False):
        lastFile = self.outputDir + sampleName + '_%s%d%d_dna_family' % (type, expectLength[0], expectLength[1])
        lastFile += '_consensus.txt'
        if (exists(lastFile)):
            print("Sequences were already analyzed " + lastFile)
            ighvSignals = {}
            faultySeqFile = self.outputDir + sampleName + '_%s%d%d_faulty_trans.fasta' % (type, expectLength[0], expectLength[1])
            noStartCodonFile = self.outputDir + sampleName + '_%s%d%d_no_atg.fasta' % (type, expectLength[0], expectLength[1])
        else:
            print("Sequences are being analyzed ... ")
            (ighvSignals, faultySeqFile, noStartCodonFile) = self.loadValidSequences(
                                            sampleName, expectLength, startCodon, type)
            
        # extract DNA motifs for each germline variant
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                self.outputDir + sampleName + 
                '_%s%d%d_dna_variant' % (type, expectLength[0], expectLength[1]),
                clusterMotifs=clusterMotifs)   
        # extract protein motifs for each each germline variant
        if expectLength[0] == expectLength[1] and type == 'secsig' :
            faultySeq = loadIGVSeqsFromFasta(faultySeqFile)
            generateMotifs(faultySeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_faulty_variant' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs) 
            noStartCodonSeq = loadIGVSeqsFromFasta(noStartCodonFile)
            generateMotifs(noStartCodonSeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_untranslated_variant' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs) 
            generateMotifs(ighvSignals, False,
                                self.outputDir + sampleName + '_%s%d%d_protein_variant' % (type, expectLength[0], expectLength[1]),
                                 transSeq=True,
                clusterMotifs=clusterMotifs)                
            
        # extract motifs for germline genes
        ighvSignals = compressSeqGeneLevel(ighvSignals)
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                            self.outputDir + sampleName + '_%s%d%d_dna_gene' % (type, expectLength[0], expectLength[1]),
                clusterMotifs=clusterMotifs) 
        if expectLength[0] == expectLength[1] and type == 'secsig':
            faultySeq = compressSeqGeneLevel(faultySeq)
            generateMotifs(faultySeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_faulty_gene' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs) 
            noStartCodonSeq = compressSeqGeneLevel(noStartCodonSeq)
            generateMotifs(noStartCodonSeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_untranslated_gene' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs)
            generateMotifs(ighvSignals, False,
                            self.outputDir + sampleName + '_%s%d%d_protein_gene' % (type, expectLength[0], expectLength[1]),
                             transSeq=True,
                clusterMotifs=clusterMotifs) 
        
        # extract motifs for germline families
        ighvSignals = compressSeqFamilyLevel(ighvSignals)
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                            self.outputDir + sampleName + '_%s%d%d_dna_family' % (type, expectLength[0], expectLength[1]),
                clusterMotifs=clusterMotifs) 
        if expectLength[0] == expectLength[1] and type == 'secsig':
            faultySeq = compressSeqFamilyLevel(faultySeq)
            generateMotifs(faultySeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_faulty_family' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs) 
            noStartCodonSeq = compressSeqFamilyLevel(noStartCodonSeq)
            generateMotifs(noStartCodonSeq, True,
                    self.outputDir + sampleName + 
                    '_%s%d%d_untranslated_family' % (type, expectLength[0], expectLength[1]),
                                 transSeq=False, extendAlphabet=True,
                clusterMotifs=clusterMotifs)
            generateMotifs(ighvSignals, False,
                            self.outputDir + sampleName + '_%s%d%d_protein_family' % (type, expectLength[0], expectLength[1]),
                             transSeq=True,
                clusterMotifs=clusterMotifs)       
        
    def refineInFramePrediction(self):
        print("In-frame prediction is being refined ...")
        queryIds = self.cloneAnnot.index
        updated = []
        procSeqs = 0
        for queryId in queryIds:
            qsRec = self.cloneAnnot.loc[queryId].to_dict()
            if qsRec['v-jframe'] == 'Out-of-frame' : # or qsRec['stopcodon'] == 'Yes'
                # unproductive sequences are excluded
                continue
#             if (queryId != '@MISEQ:151:000000000-AG0MR:1:2110:4563:8811'):
#                 continue
            inframe = True    
            # check the the v-jframe value is not NA       
            if (qsRec['v-jframe'] == 'N/A' or
                (type(qsRec['v-jframe']) != type('str') and isnan(qsRec['v-jframe']))):
#                 print "here1"
                inframe = False
            # the query sequence is not in concordance with the start of the germline gene
            qstart = qsRec['vqstart'] - qsRec['vstart'] + 1 # 1-based 
            if (not inframe or qstart < 1 or 
                (self.actualQstart != -1 and (qstart - 1 - self.actualQstart) % 3 != 0)):
#                 print "here2"
                inframe = False                
            # if no CDR3 or FR4 ==> Out-of-frame
            if (not inframe or isnan(qsRec['fr4.start']) or 
                isnan(qsRec['fr4.end'])  or isnan(qsRec['cdr3.start']) or 
               qsRec['cdr3.start'] >= qsRec['cdr3.end']):
#                 print "here3"
                inframe = False
            # doesn't start/end properly .. not multiple of 3
            if (not inframe or (qsRec['fr4.end'] - qstart + 1) % 3 != 0 or 
                (self.actualQstart != -1 and 
                 ((qsRec['fr4.end'] - self.actualQstart) % 3 != 0))):
#                 print "here4"
                inframe = False
            # indels (gaps) in FRs or CDRs cause frame-shift ==> out-of-frame
            if (not inframe or 
                qsRec['fr1.gaps'] % 3 != 0 or 
                qsRec['fr2.gaps'] % 3 != 0 or
                qsRec['fr3.gaps'] % 3 != 0 or 
                qsRec['cdr1.gaps'] % 3 != 0 or 
                qsRec['cdr2.gaps'] % 3 != 0 or
                qsRec['cdr3.gaps'] % 3 != 0):
#                 print "here5"
                inframe = False            
            # FR1 start is not aligned with query start 
#             if (not inframe or qsRec['vqstart'] != qsRec['fr1.start']):
#                 inframe = False
            if not inframe:
                self.cloneAnnot.set_value(queryId, 'v-jframe', 'Out-of-frame')
                updated += [queryId]
#             else:
#                 self.cloneAnnot.set_value(queryId, 'v-jframe', 'In-frame')
            procSeqs += 1
            if procSeqs % self.seqsPerFile == 0:
                print('%d/%d sequences have been processed ...  updated: %d' % (procSeqs, len(queryIds), len(updated)))
                sys.stdout.flush()
         
        if (len(updated) > 0):
            print("The v-j frame rearrangement has been updated for %d sequences " % len(updated))
            examples = random.choice(range(len(updated)), 10)
            for i in examples:
                print(updated[i])
            updated = None 
        # export the CDR/FR annotation to a file
        self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
        print("Text file has been written to " + self.cloneAnnotFile) 
        sys.stdout.flush()
        gc.collect()   
    
    def refineCDRInfoAnnotation(self):
        print("CDR3 annotation is being refined ...")
        # loading the 5` and 3` primers and calculate maximum alignment scores
        if self.end5:
            end5Seqs = [(rec.id, str(rec.seq), len(rec.seq)) for rec in SeqIO.parse(self.end5, "fasta")]
            L5 = max(map(lambda x:x[2], end5Seqs))
            ids = map(lambda x: x[0], end5Seqs)
            end5Seqs = map(lambda x: x[1].upper()  , end5Seqs)
            maxScores = calMaxIUPACAlignScores(end5Seqs)
            end5Seqs = zip(ids, end5Seqs, maxScores)
            valid5End = {}
            primer5End = {}
            indel5End = {}
        if self.end3:
            end3Seqs = [(rec.id, str(rec.seq), len(rec.seq)) for rec in SeqIO.parse(self.end3, "fasta")]
            L3 = max(map(lambda x:x[2], end3Seqs))
            ids = map(lambda x: x[0], end3Seqs)
            end3Seqs = map(lambda x: x[1].upper()  , end3Seqs)
            maxScores = calMaxIUPACAlignScores(end3Seqs)
            end3Seqs = zip(ids, end3Seqs, maxScores)
            valid3End = {}
            primer3End = {}
            indel3End = {}
        queryIds = self.cloneAnnot.index
        transSeqs = []
        fr1NotAtBegin = []
        endsWithStopCodon = []
        fr4NotAsExpected = []
        updatedStopCodon = 0
        stopCodonPos = {}
        noFR4 = []
        procSeqs = 0
        protein = ''
        vh = ''
        sys.stdout.flush()
        # process sequences from the FASTA file
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        else:
            records = SeqIO.index(self.readFile1, self.format)
        for id in queryIds:            
            record = records[id]
            try:    
                # retrive the sequence record from the CDRInfo file
                qsRec = self.cloneAnnot.loc[record.id].to_dict()
                seqs = [record.id, qsRec['vgene']]                
#                 if qstart <= self.actualQstart:
#                     continue
                if (qsRec['strand'] == "reversed"):
                    record.seq = record.seq.reverse_complement()
                # grab the beginning of the VH sequence  
                if self.actualQstart > -1:
                    qstart = self.actualQstart # zero-based
                else:                                  
                    qstart = int(qsRec['vqstart'] - qsRec['vstart'])  # zero-based
                if  qstart < 0:
                    qstart = 0                   
                vh = record.seq[qstart:]                       
                # check whether the VH sequence can be translated successfully
                if len(vh) % 3 != 0:
                    vh = vh[:-1 * (len(vh) % 3)]#                   
                protein = str(vh.translate())
                
                if qsRec['vqstart'] != qsRec['fr1.start']:
                    fr1NotAtBegin += [record.id]

                # FR1             
                seqs.append(extractProteinFrag(protein, qstart,
                                               qsRec['fr1.end'], qstart))
                # CDR1
                seqs.append(extractProteinFrag(protein, qsRec['cdr1.start'],
                                               qsRec['cdr1.end'], qstart))
                # FR2
                seqs.append(extractProteinFrag(protein, qsRec['fr2.start'],
                                               qsRec['fr2.end'], qstart))
                # CDR2
                seqs.append(extractProteinFrag(protein, qsRec['cdr2.start'],
                                               qsRec['cdr2.end'], qstart))
                # FR3
                seqs.append(extractProteinFrag(protein, qsRec['fr3.start'],
                                               qsRec['fr3.end'], qstart))
                # Identification of FR4 so that CDR3 can be defined 
                if isnan(qsRec['fr4.end']):
                    fr4start, fr4end = findBestAlignment(
                                extractProteinFrag(protein, qsRec['fr3.end'] + 1,
                             - 1, qstart, trimAtStop=False), FR4_CONSENSUS)#                     
                    if (fr4start != -1 and fr4end != -1 and fr4end > fr4start):
                        qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1
                        if not self.fr4cut: 
                            qsRec['fr4.end'] = len(record.seq)  # fr4end * 3 + qsRec['fr3.end']     
                        else:
                            qsRec['fr4.end'] = fr4end * 3 + qsRec['fr3.end']                                   
                        # CDR3
                        qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
                        qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
                    else:
                        qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
                        qsRec['cdr3.end'] = qsRec['jqend']
                    
                seqs.append(extractProteinFrag(protein, qsRec['cdr3.start'],
                                                   qsRec['cdr3.end'], qstart))
                seqs.append(extractProteinFrag(protein, qsRec['fr4.start'],
                                               qsRec['fr4.end'], qstart))
                ## Check whether to cut the Ig sequence after FR4 or not
                if self.fr4cut:
                    try:                        
                        protein = ''.join(seqs[2:])
                    except:
                        pass
                    try:
                        vh = record.seq[qstart:int(qsRec['fr4.end'])]
                    except:
                        pass
                if ('*' in protein):
#                     print(protein)
                    endsWithStopCodon += [record.id] 
                    # check the location of the stop codon
                    # (5-end primer, in the middle, 3-end primer)
                    stopCodonPos[record.id] = []
                    if '*' in protein[:6] :
                        stopCodonPos[record.id].append("Yes")
                    else:
                        stopCodonPos[record.id].append("No")
                    if '*' in protein[-6:]:
                        stopCodonPos[record.id].append("Yes")
                    else:
                        stopCodonPos[record.id].append("No")
                    if '*'in protein[6:-6]:
                        stopCodonPos[record.id].append("Yes")
                    else:
                        stopCodonPos[record.id].append("No")
                    ### update the StopCodon value if it was set to No
                    if qsRec['stopcodon'] == 'No':
                        updatedStopCodon += 1
                        self.cloneAnnot.set_value(record.id, 'stopcodon', 'Yes') 
                else:
                    stopCodonPos[record.id] = ["No", "No", "No"]
                # check if the primer sequences match the 5`-end and 3`-end
                
                if self.end5 and qsRec.get('5end', None) is None:
                    if self.end5offset == 0:
                        prim = str(vh[:L5])
                    else:
                        if qstart + self.end5offset >= 0:
                            prim = str(record.seq[qstart + self.end5offset: qstart + L5 + self.end5offset])
                        else:
                            prim = str(record.seq[: qstart + L5 + self.end5offset])
                    (id, tag, indelPos) = findBestMatchedPattern(prim, end5Seqs)
                    valid5End[record.id] = tag
                    primer5End[record.id] = id                    
                    indel5End[record.id] = indelPos
                if self.end3 and qsRec.get('3end', None) is None:
                    (id, tag, indelPos) = findBestMatchedPattern(str(vh[-1*L3:]), end3Seqs) 
                    valid3End[record.id] = tag
                    primer3End[record.id] = id                    
                    indel3End[record.id] = indelPos
                if (seqs[-1] != FR4_CONSENSUS):
                    fr4NotAsExpected += [record.id]
                if (seqs[-1] is None):
                    noFR4 += [record.id]
                transSeqs.append(seqs)
                #TODO: update the annotation fields with the new calculated values
                self.cloneAnnot.set_value(record.id, 'fr1.start', qstart+1)                
                gaps = abs(qsRec['vqstart'] - qsRec['vstart']) - qstart                
                mismatches = qsRec['vstart'] - 1
                if (qsRec['vstart'] > qsRec['vqstart']):
                    mismatches -= gaps
                # Only update gaps if the actual query start position is known 
                if gaps > 0:     
                    self.cloneAnnot.set_value(record.id, 'fr1.gaps', qsRec['fr1.gaps'] + gaps)                    
                    self.cloneAnnot.set_value(record.id, 'vgaps', qsRec['vgaps'] + gaps)
                # if igblast ignores mismatches at the begining ==> update
                if (mismatches > 0):
                    self.cloneAnnot.set_value(record.id, 'fr1.mismatches', qsRec['fr1.mismatches']  + mismatches)
                    self.cloneAnnot.set_value(record.id, 'vmismatches', qsRec['vmismatches']  + mismatches)
                    self.cloneAnnot.set_value(record.id, 'vstart', qsRec['vstart']  - mismatches)
                    self.cloneAnnot.set_value(record.id, 'vqstart', qsRec['vqstart']  - mismatches)
                self.cloneAnnot.set_value(record.id, 'cdr3.start', qsRec['cdr3.start'])
                self.cloneAnnot.set_value(record.id, 'cdr3.end', qsRec['cdr3.end'])
                self.cloneAnnot.set_value(record.id, 'fr4.start', qsRec['fr4.start'])
                self.cloneAnnot.set_value(record.id, 'fr4.end' , qsRec['fr4.end'])
                procSeqs += 1
                if procSeqs % self.seqsPerFile == 0:
                    print('%d/%d sequences have been processed ... ' % (procSeqs, len(queryIds)))
                    sys.stdout.flush()
            except Exception as e:                
                print("ERROR: exception in the CDR Annotation Refinement")
                print(protein, record.id, str(vh), qsRec)
                traceback.print_exc(file=sys.stdout)
                raise e

        print('%d/%d sequences have been processed ... ' % (procSeqs, len(queryIds)))
        # Expand the CDRInfo dataframe and include the 5end and 3end annotations 
        if self.end5:     
            if '5end' not in self.cloneAnnot.columns:
                self.cloneAnnot.loc[:, '5end'] = Series(valid5End, index=valid5End.keys())
                self.cloneAnnot.loc[:, '5endPrimer'] = Series(primer5End, index=primer5End.keys())
                self.cloneAnnot.loc[:, '5endIndel'] = Series(indel5End, index = indel5End.keys())

        if self.end3:    
            if '3end' not in self.cloneAnnot.columns:        
                self.cloneAnnot.loc[:, '3end'] = Series(valid3End, index=valid3End.keys())
                self.cloneAnnot.loc[:, '3endPrimer'] = Series(primer3End, index=primer3End.keys())
                self.cloneAnnot.loc[:, '3endIndel'] = Series(indel3End, index = indel3End.keys())

        ## add columns of the stop codon location 
        if 'stopat5end' not in self.cloneAnnot.columns:
            df1 = DataFrame.from_dict(stopCodonPos, orient='index')
            df1.columns = ['stopat5end', 'stopat3end', 'stopinmiddle']
            self.cloneAnnot = concat([self.cloneAnnot, df1], axis=1)

        ## print statsitics and final processed data 
        if (len(fr1NotAtBegin) > 0):
            print("%d sequences have FR1 start not equal to query start (Excluded)" % (len(fr1NotAtBegin)))
            examples = random.choice(range(len(fr1NotAtBegin)), 10)
            for i in examples:
                print(fr1NotAtBegin[i])
            fr1NotAtBegin = None 
        if (len(endsWithStopCodon) > 0):
            print("%d sequences contain a stop codon " % (len(endsWithStopCodon)))
            examples = random.choice(range(len(endsWithStopCodon)), 10)
            for i in examples:
                print(endsWithStopCodon[i])
            endsWithStopCodon = None   
        if (updatedStopCodon > 0):
            print("The stopcodon flag was updated for %d sequencs " % (updatedStopCodon))          
        if (len(fr4NotAsExpected) > 0):
            print("%d sequences do not have an expected FR4 sequences (%s) " % (len(fr4NotAsExpected), FR4_CONSENSUS))
            examples = random.choice(range(len(fr4NotAsExpected)), 10)
            for i in examples:
                print(fr4NotAsExpected[i])
            fr4NotAsExpected = None
        if (len(noFR4) > 0):
            print("%d sequences do not have FR4 " % (len(noFR4)))
            examples = random.choice(range(len(noFR4)), 10)
            for i in examples:
                print(noFR4[i])
            noFR4 = None
            
        self.cdrSeqs = DataFrame(transSeqs, columns=['queryid', 'germline', 'fr1', 'cdr1', 'fr2', 'cdr2',
                                                  'fr3', 'cdr3', 'fr4'])
        self.cdrSeqs.index = self.cdrSeqs.queryid
        del self.cdrSeqs['queryid']
        self.cdrSeqs.to_csv(self.cdrSeqFile, sep='\t', header=True, index=True)
        
        # export the CDR/FR annotation to a file
        self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
        print("Text file has been written to " + self.cloneAnnotFile) 
        sys.stdout.flush()
        gc.collect()
     
    def analyzeRestrictionSitesSimple(self):        
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        outDir = self.outputDir + "restriction_sites/"
        if (not os.path.isdir(outDir)):
            os.system("mkdir " + outDir)
        self.siteHitsFile = outDir + self.name
        self.siteHitsFile += "_%s_simple.csv" % (self.sitesFile.split('/')[-1].split('.')[0]) 
        
        if (exists(self.siteHitsFile)):
            print("Restriction sites were already searched at ... " + self.siteHitsFile.split('/')[-1])
            return      
               
        self.loadRestrictionSites()
        print("Restriction sites are being searched ... ") 
        sys.stdout.flush()
        siteHitsCount = {}
        siteHitSeqsCount = {}        
        siteHitSeqsGermline = {}
        seqsCutByAny = 0
        siteHitsSeqsIDs = {}
#         siteHitsSeqsIGV = {}
        for site in self.sites.keys():
            siteHitsCount[site] = 0
            siteHitSeqsCount[site] = 0            
            siteHitSeqsGermline[site] = []
            siteHitsSeqsIDs[site] = set()
#             siteHitsSeqsIGV[site] = set()
        
        procSeqs = 0
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        else:
            records = SeqIO.index(self.readFile1, self.format)
        for id in records.keys():
            record = records[id]
            try:              
                seq = str(record.seq)
                seqRC = str(Seq(seq).reverse_complement())
                cut = False
                for site in siteHitsCount.keys():
                    hits = findHits(seq, self.sites[site])
                    strand = "forward"
                    if len(hits) == 0:
                        hits = findHits(seqRC, self.sites[site])
                        strand = "reversed"
                    if len(hits) > 0:
                        siteHitsCount[site] += len(hits) 
                        siteHitSeqsCount[site] += 1                     
                        siteHitsSeqsIDs[site].add(record.id)   
                        cut = True                 
                if cut:
                    seqsCutByAny += 1
                procSeqs += 1
                if procSeqs % self.seqsPerFile == 0:
                    print('%d/%d sequences have been searched ... ' % (procSeqs, len(records.keys())))
                    sys.stdout.flush()
#                 break
            except BaseException as e:                
                print(e)
                raise
        print('%d/%d sequences have been searched ... ' % (procSeqs, len(records.keys())))
        # # print out the results
        f = open(self.siteHitsFile, 'w')
        f.write("Enzyme,Restriction Site,No.Hits,Percentage of Hits (%),No.Molecules,Percentage of Molecules (%) \n")
        sites = sorted(siteHitSeqsCount, key=siteHitSeqsCount.get)
        for site in sites:
            f.write("%s,%s,%d,%.3f,%d,%.3f\n" % (site ,
                                              self.sites[site],
                                              siteHitsCount[site],
                                            siteHitsCount[site] * 100.0 / sum(siteHitsCount.values()),
                                            siteHitSeqsCount[site],
                                            siteHitSeqsCount[site] * 100.0 / len(records.keys())))
            # write the first 100 sequences cut in the germline of each restriction enzyme 
            
        f.write("Sequences cut by any of the above enzymes, , , , %d, %.3f\n" % (seqsCutByAny, seqsCutByAny * 100.0 / len(records.keys())))
        f.close()
        # Ven Diagram of overlapping sequences
        plotVenn(siteHitsSeqsIDs, self.siteHitsFile.replace('.csv', '_venn.png'))
        print("Restriction enzyme results were written to " + self.siteHitsFile)
          
       
    def analyzeRestrictionSites(self):        
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        self.siteHitsFile = self.outputDir + self.name
        self.siteHitsFile += "_%s.csv" % (self.sitesFile.split('/')[-1].split('.')[0]) 
        
        if (exists(self.siteHitsFile)):
            print("Restriction sites were already searched at ... " + self.siteHitsFile.split('/')[-1])
            return
        
        self.cloneAnnotFile = self.outputDir + self.name
        self.cloneAnnotFile += "_cdr_info.tab" 
        
        self.cdrSeqFile = self.outputDir + self.name
        self.cdrSeqFile += "_cdr_seq.tab"
        if (not exists(self.cloneAnnotFile)):       
            self.analyzeAbundance('cdrinfo')  
            # export the CDR/FR annotation to a file
            self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
            print("Text file has been written to " + self.cloneAnnotFile) 
            # Refine the annotation of CDR3 and FR4
            self.refineCDRInfoAnnotation()  
        else:
            print("\tCDRs information file was found! ... " + self.cloneAnnotFile.split('/')[-1])
            self.cloneAnnot = read_csv(self.cloneAnnotFile, sep='\t',
                                       header=0, index_col=0)            
            if (not exists(self.cdrSeqFile)):                
                self.refineCDRInfoAnnotation()  
        
        self.loadRestrictionSites()
        print("Restriction sites are being searched ... ")
        self.cdrSeqs = None
        gc.collect()
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'No']
        queryIds = self.cloneAnnot.index
        siteHitsCount = {}
        siteHitSeqsCount = {}
        hitRegion = {}
        siteHitSeqsGermline = {}
        seqsCutByAny = 0
        siteHitsSeqsIDs = {}
        siteHitsSeqsIGV = {}
        for site in self.sites.keys():
            siteHitsCount[site] = 0
            siteHitSeqsCount[site] = 0
            hitRegion[site] = Counter({'fr1':0, 'cdr1':0,
                                       'fr2':0, 'cdr2':0,
                                       'fr3':0, 'cdr3':0,
                                       'fr4':0})
            siteHitSeqsGermline[site] = []
            siteHitsSeqsIDs[site] = set()
            siteHitsSeqsIGV[site] = set()
        germline = set(['fr1', 'fr2', 'fr3', 'cdr1', 'cdr2'])
        procSeqs = 0
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        else:
            records = SeqIO.index(self.readFile1, self.format)
        for id in queryIds:
            record = records[id]
            try:
                qsRec = self.cloneAnnot.loc[record.id].to_dict()
                qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
                if (isnan(qsRec['fr4.end'])):
                    end = len(record.seq)
                else:
                    end = int(qsRec['fr4.end'])
                seq = str(record.seq[qstart:end])
                seqRC = str(Seq(seq).reverse_complement())
                cut = False
                for site in siteHitsCount.keys():
                    hits = findHits(seq, self.sites[site])
                    strand = "forward"
                    if len(hits) == 0:
                        hits = findHits(seqRC, self.sites[site])
                        strand = "reversed"
                    if len(hits) > 0:
                        siteHitsCount[site] += len(hits) 
                        siteHitSeqsCount[site] += 1
                        hitsRegion = findHitsRegion(qsRec, hits)
                        if (len(set(hitsRegion).intersection(germline)) > 0 
                            and len(siteHitSeqsGermline[site]) < 10000):
                            siteHitSeqsGermline[site].append((strand, str(record.seq)))  
                            siteHitsSeqsIGV[site].add(qsRec['vgene'].split('*')[0])                          
                        hitRegion[site] += Counter(hitsRegion)
                        siteHitsSeqsIDs[site].add(record.id)   
                        cut = True                 
                if cut:
                    seqsCutByAny += 1
                procSeqs += 1
                if procSeqs % self.seqsPerFile == 0:
                    print('%d/%d sequences have been searched ... ' % (procSeqs, len(queryIds)))
#                 break
            except BaseException as e:
                print(qstart, end, len(record.seq), str(record.seq))
                print(e)
                raise
        print('%d/%d sequences have been searched ... ' % (procSeqs, len(queryIds)))
        # # print out the results
        f = open(self.siteHitsFile, 'w')
        f.write("Enzyme,Restriction Site,No.Hits,Percentage of Hits (%),No.Molecules,Percentage of Molecules (%),FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4, V Germlines \n")
        sites = sorted(siteHitSeqsCount, key=siteHitSeqsCount.get)
        for site in sites:
            f.write("%s,%s,%d,%.3f,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%s\n" % (site ,
                                              self.sites[site],
                                              siteHitsCount[site],
                                            siteHitsCount[site] * 100.0 / sum(siteHitsCount.values()),
                                            siteHitSeqsCount[site],
                                            siteHitSeqsCount[site] * 100.0 / len(queryIds),
                                            hitRegion[site]['fr1'],
                                            hitRegion[site]['cdr1'],
                                            hitRegion[site]['fr2'],
                                            hitRegion[site]['cdr2'],
                                            hitRegion[site]['fr3'],
                                            hitRegion[site]['cdr3'],
                                            hitRegion[site]['fr4'],
                                            '|'.join(siteHitsSeqsIGV[site])))
            # write the first 100 sequences cut in the germline of each restriction enzyme 
            seqs = []
            for (strand, seq) in siteHitSeqsGermline[site]:
                seqs.append(SeqRecord(Seq(seq), id='seq' + `len(seqs)` + strand))
            SeqIO.write(seqs, self.siteHitsFile.replace('.csv', '_germline' + site + '.fasta'), 'fasta')
        f.write("Sequences cut by any of the above enzymes, %d, %.3f\n" % (seqsCutByAny, seqsCutByAny * 100.0 / len(queryIds)))
        f.close()
        # Ven Diagram of overlapping sequences
        plotVenn(siteHitsSeqsIDs, self.siteHitsFile.replace('.csv', '_venn.png'))
        print("Restriction enzyme results were written to " + self.siteHitsFile)
        
           
    def loadRestrictionSites(self):
        f = open(self.sitesFile)
        lines = f.readlines()
        self.sites = {}
        for line in lines:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    fields = line.split()
                    if (self.sites.get(fields[0], None) is not None):
                        print(fields[0] + " is duplicated.")
                    site = fields[1].upper()
                    site = replaceIUPACLetters(site)
                    site = site.replace('N', '.').replace('(', '[').replace(')', ']')                
                    self.sites[fields[0]] = site
                except:
                    print line, fields
                    raise
        print("\t\tRestricting sites have been loaded")
                
    def analyzeDiversity(self):
        print("The diversity of the CDRs of the VH regions is being analyzed ... ")
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]      
        if (not os.path.isdir(self.outputDir + "cdrs/")):
            os.system("mkdir " + self.outputDir + "cdrs/")
        self.cloneAnnotFile = self.outputDir + "cdrs/" + self.name
        self.cloneAnnotFile += "_cdr_info.tab"   
        
        self.cdrSeqFile = self.outputDir + "cdrs/" + self.name
        self.cdrSeqFile += "_cdr_seq.tab"     
        
        if (not exists(self.cloneAnnotFile)):       
            self.analyzeAbundance('cdrinfo')
            # export the CDR/FR annotation to a file
            self.cloneAnnot.to_csv(self.cloneAnnotFile, sep='\t', header=True, index=True)    
            print("Text file has been written to " + self.cloneAnnotFile) 
            # Refine the annotation of CDR3 and FR4
            self.refineCDRInfoAnnotation()  
            self.refineInFramePrediction()
            
        else:
            print("\tCDRs information file was found! ... " + self.cloneAnnotFile.split('/')[-1])
            sys.stdout.flush()
            self.cloneAnnot = read_csv(self.cloneAnnotFile, sep='\t',
                                       header=0, index_col=0)            
            if (not exists(self.cdrSeqFile)):                
                self.refineCDRInfoAnnotation() 
                self.refineInFramePrediction() 
            else:
                print("\tCDRs sequence file was found! ... " + self.cdrSeqFile.split('/')[-1])
                sys.stdout.flush()
                self.cdrSeqs = read_csv(self.cdrSeqFile, sep='\t',
                                       header=0, index_col=0)            
            
                 
#         self.alignInfo = read_csv(self.cloneAnnotFile.replace('cdr', 'align'), 
#                                   sep ='\t', header= 0, index_col=0)
#         cdrids = set(self.cloneAnnot.index)
#         alignids = set(self.alignInfo.index)
#         incdr = cdrids.difference(alignids)
#         inalign = alignids.difference(cdrids) 
#         print(len(cdrids), len(alignids))
#         print(len(incdr), len(inalign))
#         print(list(inalign)[1:10])
        
        outDir = self.outputDir
        self.outputDir = outDir + "/cdrs/"
                
        self.extractProductiveRNAs()
        sys.stdout.flush()
        
        self.cdrSeqs = self.cdrSeqs[map(lambda x: x in self.cloneAnnot.index, self.cdrSeqs.index)]
                
        self.writeCDRStats()
        sys.stdout.flush()
        
        self.writeFRStats()
        sys.stdout.flush()
        
        self.writeVGeneStats()
        sys.stdout.flush()
        
        self.generateCDRandFRLogos()
        sys.stdout.flush()
        # analyze productive clones only 
        self.analyzeIgProtein()
        sys.stdout.flush()
        
        self.outputDir = outDir
        
    def analyzeIgProtein(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
#         self.cdrSeqs = read_csv(self.cdrSeqFile, sep='\t',
#                                        header=0, index_col=0)
        self.readFile1 = self.outputDir + self.name
        self.readFile1 += '_productive_prot.fasta'
        if (not exists(self.readFile1)):
            print("Protein sequences are being prepared ...")
            records = []
            procSeqs = 0
            open(self.readFile1, 'w').close()
            for id in self.cdrSeqs.index:
                seq = ''.join(self.cdrSeqs.loc[id, ].tolist()[1:])
                if '*' in seq:
                    seq = seq.replace('*', 'X')
                rec = SeqRecord(Seq(seq), id=id, name="", description="")
                records.append(rec)
                procSeqs += 1              
                if procSeqs % self.seqsPerFile == 0:
                    print('\t%d/%d sequences have been processed ...  ' % (procSeqs, len(self.cdrSeqs.index)))
                    sys.stdout.flush()
                    SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
                    records = []    
            SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
            del records
        else:
            print("File found ... " + self.readFile1.split('/')[-1])
        self.format = 'fasta'
        self.readFile2 = None
        self.seqType = 'protein'
        self.bitScore = [0, Inf] 
        self.alignLen = [0, Inf]
        self.sStart = [1, Inf]
        if (exists(self.outputDir + "/abundance/")):
            print("Protein sequences have been already analyzed ... ")
        else:
            self.analyzeAbundance('abundance')
    
    def write5EndPrimerStats(self, cdrInfo, fileprefix, category="All"):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        
        valid5End = Counter(cdrInfo['5end'].tolist())
        plotDist(valid5End, self.name, fileprefix + 'integrity_dist.png', 
                 title='Integrity of 5`-end Primer Sequence (%s)' % (category),
             proportion=True, rotateLabels=False)     
        invalid5Clones = cdrInfo.index[cdrInfo['5end'] == 'Indelled'].tolist()
        print("Example of Indelled 5`-end:", invalid5Clones[1:10])
        print("Example of valid 5`-end:", cdrInfo.index[cdrInfo['5end'] != 'Indelled'].tolist()[1:10])
        
        stopcodonInFrameDist = Counter(cdrInfo['stopcodon'].tolist())
        plotDist(stopcodonInFrameDist, self.name, 
                 fileprefix + 'stopcodon_dist.png', 
                 title='Stop Codons in 5`-End (%s)' % (category),
                 proportion=False, rotateLabels=False)
        
        c1 = Counter(cdrInfo[cdrInfo['5end'] == 'Indelled']['5endPrimer'].tolist())
        plotDist(c1, self.name, fileprefix + 
             'indelled_dist.png', 
             title='Abundance of Indelled 5`-end Primers (%s)' % (category),
             proportion=False, rotateLabels=False, vertical=False, top=50)
        c = Counter(cdrInfo[cdrInfo['5end'] == 'Indelled']['5endIndel'].tolist())
        plotDist(c, self.name, fileprefix + 
             'indel_pos_dist.png', 
             title='Abundance of Indel Positions in 5`-end Primers (%s)' % (category),
             proportion=False, rotateLabels=False, vertical=True, 
             sortValues=False, top=50)
        primers = set(cdrInfo['5endPrimer'].tolist())
#         print(c1, primers)
        for primer in primers:
#             print(primer)
            df = cdrInfo[cdrInfo['5end'] == 'Indelled']
            df = df[df['5endPrimer'] == primer]
#             print(df.shape)
            germLineDist = compressCountsGeneLevel(Counter(df['vgene'].tolist()))        
            plotDist(germLineDist, self.name, fileprefix + primer +  
                     '_igv_dist.png',
                      title='IGV Abundance (%s)' % (category),
                     proportion=False, vertical=False, top=20, rotateLabels=False)
        gc.collect()
        
    def generateCDRandFRLogos(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        seqs = {}
        for k in self.cdrSeqs.columns:
            if k.startswith('cdr') :
                seqs[k] = self.cdrSeqs[k].tolist()
            if k.startswith('fr'):
                seqs[k] = self.cdrSeqs[k].tolist()
        # generate Toby's logos ?!?!
        generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_fr_Toby')
        
        # generate CDR/FR logos without alignment
        generateMotifs(seqs, False, self.outputDir + self.name + '_cdr_fr',
                       protein=True) 
        # generate CDR/FR logos after alignment
        generateMotifs(seqs, True,
                    self.outputDir + self.name + '_cdr_fr_aligned',
                    protein=True)
        # generate CDR3 logos per  germline 
        seqs = {}
        vgenes = self.cdrSeqs["germline"].tolist()
        vgenes = map(lambda x: x.split('*')[0], vgenes)
        for vgene in set(vgenes):
            seqs["cdr3_"+vgene.replace("/", "-")] = self.cdrSeqs.loc[map(lambda x: x == vgene, vgenes), "cdr3"].tolist()            
        # generate Toby's logos ?!?!
        generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_Toby_gene')
        # generate CDR/FR logos after alignment
        generateMotifs(seqs, True,
                    self.outputDir + self.name + '_cdr_gene_aligned',
                    protein=True)
        # generate CDR3 logos per family        
        seqs = {}
        vfams = map(lambda x: x.split('-')[0].split('/')[0], vgenes)
        for vfam in set(vfams):
            seqs["cdr3_"+vfam] = self.cdrSeqs.loc[map(lambda x: x == vgene, vgenes), "cdr3"].tolist()
        # generate Toby's logos ?!?!
        generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_fr_Toby')
        # generate CDR/FR logos after alignment
        generateMotifs(seqs, True,
                    self.outputDir + self.name + '_cdr_fr_aligned',
                    protein=True)
        
       
    def writeVGeneStats(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        # CDR1 statistics 
        vLength = (self.cloneAnnot['fr3.end'] - self.cloneAnnot['fr1.start'] + 1) / 3
        vLength = vLength.tolist()
        histcals = plotSeqLenDist(vLength, self.name, self.outputDir + self.name + 
                 '_vgene_len_dist.png', dna=False,
                  seqName='V Gene', normed=True, maxbins=20)
#         print(histcals)
        vGaps = Counter(self.cloneAnnot['vgaps'].tolist())
        plotDist(vGaps, self.name, self.outputDir + self.name + 
                 '_vgene_gaps_dist.png', title='Gaps in V Gene',
                 proportion=True, rotateLabels=False, top=20) 
        cdrMismatches = Counter(self.cloneAnnot['vmismatches'].tolist())
        plotDist(cdrMismatches, self.name, self.outputDir + self.name + 
                 '_vgene_mismatches_dist.png', title='Mismatches in V Gene',
                 proportion=True, rotateLabels=False, top=20) 
        # quantify V domain sequence diversity                
        if (not exists(self.outputDir + self.name + 
                 '_Vdomain_duplication_family.png')):
            print("Grouping V domain sequences per family ...")
            VH = {}
    #         i = 0
            ighvs = map(lambda x : x.split('-')[0].split('/')[0], self.cdrSeqs['germline'].tolist())
            for ighv in set(ighvs):
                VH[ighv] = []
            for (ighv, f1, c1, f2, c2, f3, c3, f4) in zip(ighvs,
                                                        self.cdrSeqs['fr1'].tolist(),
                                                          self.cdrSeqs['cdr1'].tolist(),
                                                          self.cdrSeqs['fr2'].tolist(),
                                                          self.cdrSeqs['cdr2'].tolist(),
                                                          self.cdrSeqs['fr3'].tolist(),
                                                          self.cdrSeqs['cdr3'].tolist(),
                                                          self.cdrSeqs['fr4'].tolist()):           
                try:
                    VH[ighv].append(''.join([f1, c1, f2, c2, f3, c3, f4]))
                except:
                    if (f4 is None or isnan(f4)):  # or c3 is None or isnan(c3):
                        VH += [''.join([f1, c1, f2, c2, f3, c3])]
                    else:
                        print(id, f1, c1, f2, c2, f3, c3, f4)
            ighvs = VH.keys()
            ighvs.sort()        
            
            plotSeqDuplication(map(lambda x:VH[x], ighvs),
                             self.outputDir + self.name + 
                     '_Vdomain_duplication_family.png',
                             ighvs,
                             'Duplication of V Domain Sequences Per Family', True)
            gc.collect()
        
    
    def extractProductiveRNAs(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        # v-j rearrangement frame distribution 
        vjframeDist = Counter(self.cloneAnnot['v-jframe'].tolist())
        if NaN in vjframeDist.keys():
            nanCounts = vjframeDist[NaN]
            vjframeDist = Counter({'In-frame': vjframeDist['In-frame'],
                                   'Out-of-frame': vjframeDist['Out-of-frame'] + nanCounts})
        plotDist(vjframeDist, self.name, self.outputDir + self.name + 
                 '_vjframe_dist.png', title='V-D-J Rearrangement',
                 proportion=False, rotateLabels=False)
        print(vjframeDist)
        del vjframeDist
        if self.end5:
            print("5-end analysis of all clones ... ")
            self.write5EndPrimerStats(self.cloneAnnot, self.outputDir+self.name+
                                      '_all_5end_')
            invalid5Clones = self.cloneAnnot.index[self.cloneAnnot['5end'] == 'Indelled'].tolist()
        if self.end3:
            valid3End = Counter(self.cloneAnnot['3end'].tolist())
            plotDist(valid3End, self.name, self.outputDir + self.name + 
                 '_all_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence',
                 proportion=True, rotateLabels=False)
            invalid3Clones = self.cloneAnnot.index[self.cloneAnnot['3end'] == 'Indelled'].tolist()
            print("Example of Indelled 3`-end:", invalid3Clones[1:10])
            try:
                plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
                          self.outputDir + self.name + 
                     '_all_invalid_primers.png')
                del invalid5Clones, invalid3Clones
            except:
                pass
            del valid3End
        
        OutOfFrame = self.cloneAnnot[self.cloneAnnot['v-jframe'] != 'In-frame']
        OutOfFrameFamilyDist = compressCountsFamilyLevel(Counter(OutOfFrame['vgene'].tolist()))
        plotDist(OutOfFrameFamilyDist, self.name, self.outputDir + self.name + 
                 '_notinframe_igv_dist.png',
                  title='IGV Abundance of Not In-frame Sequences',
                 proportion=True)
        del OutOfFrameFamilyDist
        OutOfFrame = OutOfFrame[OutOfFrame['v-jframe'] == 'Out-of-frame']
        cdrLength = (OutOfFrame['cdr1.end'] - OutOfFrame['cdr1.start'] + 1) / 3
        cdrLength = cdrLength.tolist()
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_outframe_cdr1_len_dist.png', dna=False,
                  seqName='CDR1', normed=True, maxbins=10)
        cdrGaps = Counter(OutOfFrame['cdr1.gaps'].tolist())
        plotDist(cdrGaps, self.name, self.outputDir + self.name + 
                 '_outframe_cdr1_gaps_dist.png', title='Gaps in CDR1',
                 proportion=False, rotateLabels=False)
        frGaps = Counter(OutOfFrame['fr1.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_outframe_fr1_gaps_dist.png', title='Gaps in FR1',
                 proportion=False, rotateLabels=False)
        del cdrLength, cdrGaps, frGaps
        cdrLength = (OutOfFrame['cdr2.end'] - OutOfFrame['cdr2.start'] + 1) / 3
        cdrLength = cdrLength.tolist()
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_outframe_cdr2_len_dist.png', dna=False,
                  seqName='CDR2', normed=True, maxbins=10)
        cdrGaps = Counter(OutOfFrame['cdr2.gaps'].tolist())
        plotDist(cdrGaps, self.name, self.outputDir + self.name + 
                 '_outframe_cdr2_gaps_dist.png', title='Gaps in CDR2',
                 proportion=False, rotateLabels=False)
        frGaps = Counter(OutOfFrame['fr2.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_outframe_fr2_gaps_dist.png', title='Gaps in FR2',
                 proportion=False, rotateLabels=False)
        del cdrLength, cdrGaps, frGaps
        cdrGaps = Counter([x if not isnan(x) else 'NA' for x in OutOfFrame['cdr3.gaps'] ])
#         print(len(cdrGaps))
        plotDist(cdrGaps, self.name, self.outputDir + self.name + 
                 '_outframe_cdr3_gaps_dist.png', title='Gaps in Germline CDR3',
                 proportion=False, rotateLabels=False)
        frGaps = Counter(OutOfFrame['fr3.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_outframe_fr3_gaps_dist.png', title='Gaps in FR3',
                 proportion=False, rotateLabels=False)
        del cdrGaps, frGaps
        if self.end5:
            print("5-end analysis of out-of-frame clones ... ")
            self.write5EndPrimerStats(OutOfFrame, self.outputDir+self.name+
                                      '_outframe_5end_', 'Out-of-frame') 
            invalid5Clones = OutOfFrame.index[OutOfFrame['5end'] == 'Indelled'].tolist()                     
        if self.end3:
            valid3End = Counter(OutOfFrame['3end'].tolist())
            plotDist(valid3End, self.name, self.outputDir + self.name + 
                 '_outframe_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence(Out-of-frame)',
                 proportion=True, rotateLabels=False)
            invalid3Clones = OutOfFrame.index[OutOfFrame['3end'] == 'Indelled'].tolist()
            print("Example of out-of-frame Indelled 3`-end:", invalid3Clones[1:10])
            print("Example of out-of-frame valid 3`-end:", OutOfFrame.index[OutOfFrame['3end'] != 'Indelled'].tolist()[1:10])
            try:
                plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
                          self.outputDir + self.name + 
                     '_outframe_invalid_primers.png')
                del invalid5Clones, invalid3Clones
            except Exception as e:
                raise e
            del valid3End
        del OutOfFrame
        # choose only In-frame RNA sequences
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']
        # Stop Codon 
        stopcodonInFrameDist = Counter(self.cloneAnnot['stopcodon'].tolist())
        plotDist(stopcodonInFrameDist, self.name, self.outputDir + self.name + 
                 '_inframe_stopcodon_dist.png', title='In-frame Stop Codons',
                 proportion=False, rotateLabels=False)
        print(stopcodonInFrameDist)
        # stop codon family distribution
        stopcodFamily = Counter(self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'Yes']['vgene'].tolist())
        stopcodFamily = compressCountsFamilyLevel(stopcodFamily)
        plotDist(stopcodFamily, self.name, self.outputDir + self.name + 
                 '_inframe_stopcodfam_dist.png',
                  title='IGV Abundance of In-frame Unproductive Sequences',
                 proportion=True)
        del stopcodonInFrameDist, stopcodFamily
#         print(stopcodFamily)
        # choose only productive RNA sequences 
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'No']
        productiveFamilyDist = compressCountsFamilyLevel(Counter(self.cloneAnnot['vgene'].tolist()))
        plotDist(productiveFamilyDist, self.name, self.outputDir + self.name + 
                 '_productive_igv_dist.png',
                  title='IGV Abundance of Productive Sequences',
                 proportion=True)
        del productiveFamilyDist
        if self.end5:
            valid5End = Counter(self.cloneAnnot['5end'].tolist())
            plotDist(valid5End, self.name, self.outputDir + self.name + 
                 '_productive_5end_integrity_dist.png', title='Integrity of 5`-end Primer Sequence(Productive)',
                 proportion=True, rotateLabels=False)
            invalid5Clones = self.cloneAnnot.index[self.cloneAnnot['5end'] == 'Indelled'].tolist()
            print("Example of invalid 5`-end:", invalid5Clones[1:10])
        if self.end3:
            valid3End = Counter(self.cloneAnnot['3end'].tolist())
            plotDist(valid3End, self.name, self.outputDir + self.name + 
                 '_productive_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence(Productive)',
                 proportion=True, rotateLabels=False)
            invalid3Clones = self.cloneAnnot.index[self.cloneAnnot['3end'] == 'Indelled'].tolist()
            print("Example of invalid 3`-end:", invalid3Clones[1:10])
            try:
                plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
                          self.outputDir + self.name + 
                     '_productive_invalid_primers.png')
            except Exception as e:
                raise e
        gc.collect()
        
    def writeFRStats(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]

        # FR1 statistics 
        frLength = (self.cloneAnnot['fr1.end'] - self.cloneAnnot['fr1.start'] + 1) / 3
        frLength = frLength.tolist()
        histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
                 '_fr1_len_dist.png', dna=False,
                  seqName='FR1', normed=True, maxbins=20)
#         print(histcals)
        frGaps = Counter(self.cloneAnnot['fr1.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_fr1_gaps_dist.png', title='Gaps in FR1',
                 proportion=False, rotateLabels=False) 
        frMismatches = Counter(self.cloneAnnot['fr1.mismatches'].tolist())
        plotDist(frMismatches, self.name, self.outputDir + self.name + 
                 '_fr1_mismatches_dist.png', title='Mismatches in FR1',
                 proportion=False, rotateLabels=False) 
        # FR2 statistics 
        frLength = (self.cloneAnnot['fr2.end'] - self.cloneAnnot['fr2.start'] + 1) / 3
        frLength = frLength.tolist()
        histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
                 '_fr2_len_dist.png', dna=False,
                  seqName='FR2', normed=True, maxbins=20)
#         print(histcals)
        frGaps = Counter(self.cloneAnnot['fr2.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_fr2_gaps_dist.png', title='Gaps in FR2',
                 proportion=False, rotateLabels=False) 
        frMismatches = Counter(self.cloneAnnot['fr2.mismatches'].tolist())
        plotDist(frMismatches, self.name, self.outputDir + self.name + 
                 '_fr2_mismatches_dist.png', title='Mismatches in FR2',
                 proportion=False, rotateLabels=False) 
        # FR3 statistics 
        frLength = (self.cloneAnnot['fr3.end'] - self.cloneAnnot['fr3.start'] + 1) / 3
        frLength = frLength.tolist()
        histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
                 '_fr3_len_dist.png', dna=False,
                  seqName='FR3', normed=True, maxbins=20)
#         print(histcals)
        frGaps = Counter(self.cloneAnnot['fr3.gaps'].tolist())
        plotDist(frGaps, self.name, self.outputDir + self.name + 
                 '_fr3_gaps_dist.png', title='Gaps in FR3',
                 proportion=False, rotateLabels=False) 
        frMismatches = Counter(self.cloneAnnot['fr3.mismatches'].tolist())
        plotDist(frMismatches, self.name, self.outputDir + self.name + 
                 '_fr3_mismatches_dist.png', title='Mismatches in FR3',
                 proportion=False, rotateLabels=False)
        # FR4
        frLength = (self.cloneAnnot['fr4.end'] - self.cloneAnnot['fr4.start'] + 1) / 3
        frLength = frLength.tolist()
        histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
                 '_fr4_len_dist.png', dna=False,
                  seqName='FR4', normed=True, maxbins=20)
        gc.collect()
        # Quantify FR sequence diversity
        plotSeqDuplication([self.cdrSeqs['fr1'].tolist(),
                          self.cdrSeqs['fr2'].tolist(),
                          self.cdrSeqs['fr3'].tolist(),
                          self.cdrSeqs['fr4'].tolist()],
                         self.outputDir + self.name + 
                 '_fr_duplication.png',
                         ['FR1', 'FR2', 'FR3', 'FR4'],
                         'Duplication of FR Sequences')
        gc.collect()
        plotSeqDiversity([self.cdrSeqs['fr1'].tolist(),
                          self.cdrSeqs['fr2'].tolist(),
                          self.cdrSeqs['fr3'].tolist(),
                          self.cdrSeqs['fr4'].tolist()],
                         self.outputDir + self.name + 
                 '_fr_diversity.png',
                         ['FR1', 'FR2', 'FR3', 'FR4'],
                         'Diversity of FR Sequences')
        
    def writeCDRStats(self):
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        # CDR1 statistics 
        cdrLength = (self.cloneAnnot['cdr1.end'] - self.cloneAnnot['cdr1.start'] + 1) / 3
        cdrLength = cdrLength.tolist()
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_cdr1_len_dist.png', dna=False,
                  seqName='CDR1', normed=True, maxbins=20)
#         print(histcals)
        cdrGaps = Counter(self.cloneAnnot['cdr1.gaps'].tolist())
        plotDist(cdrGaps, self.name, self.outputDir + self.name + 
                 '_cdr1_gaps_dist.png', title='Gaps in CDR1',
                 proportion=False, rotateLabels=False) 
        cdrMismatches = Counter(self.cloneAnnot['cdr1.mismatches'].tolist())
        plotDist(cdrMismatches, self.name, self.outputDir + self.name + 
                 '_cdr1_mismatches_dist.png', title='Mismatches in CDR1',
                 proportion=False, rotateLabels=False) 
#         print(cdrGaps)
        # CDR2 stats
        cdrLength = (self.cloneAnnot['cdr2.end'] - self.cloneAnnot['cdr2.start'] + 1) / 3
        cdrLength = cdrLength.tolist()
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_cdr2_len_dist.png', dna=False,
                 seqName='CDR2', normed=True, maxbins=20)
#         print(histcals)
        cdrGaps = Counter(self.cloneAnnot['cdr2.gaps'].tolist())
        plotDist(cdrGaps, self.name, self.outputDir + self.name + 
                 '_cdr2_gaps_dist.png', title='Gaps in CDR2',
                 proportion=False, rotateLabels=False)
#         print(cdrGaps)
        cdrMismatches = Counter(self.cloneAnnot['cdr2.mismatches'].tolist())
        plotDist(cdrMismatches, self.name, self.outputDir + self.name + 
                 '_cdr2_mismatches_dist.png', title='Mismatches in CDR2',
                 proportion=False, rotateLabels=False)
        # CDR3 stats
        cdrLength = (self.cloneAnnot['cdr3.end'] - self.cloneAnnot['cdr3.start'] + 1) / 3
        cdrLength = cdrLength.tolist()
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_cdr3_len_dist_nooutliers.png', dna=False,
                 seqName='CDR3', normed=True, maxbins=20,
                 removeOutliers=True)
        histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
                 '_cdr3_len_dist.png', dna=False,
                 seqName='CDR3', normed=True, maxbins=20,
                 removeOutliers=False)
        gc.collect()
#         print(histcals)
        # quantify CDR sequence diversity      
        
        if (not exists(self.outputDir + self.name + 
                 '_Vdomain_diversity.png')):            
    #         i = 0
            VH = []
            for (id, f1, c1, f2, c2, f3, c3, f4) in zip(self.cdrSeqs.index.tolist(),
                                                        self.cdrSeqs['fr1'].tolist(),
                                                          self.cdrSeqs['cdr1'].tolist(),
                                                          self.cdrSeqs['fr2'].tolist(),
                                                          self.cdrSeqs['cdr2'].tolist(),
                                                          self.cdrSeqs['fr3'].tolist(),
                                                          self.cdrSeqs['cdr3'].tolist(),
                                                          self.cdrSeqs['fr4'].tolist()):           
                try:
                    VH += [''.join([f1, c1, f2, c2, f3, c3, f4])]
                except:
                    if (f4 is None or isnan(f4)):  # or c3 is None or isnan(c3):
                        VH += [''.join([f1, c1, f2, c2, f3, c3])]
                    else:
                        print(id, f1, c1, f2, c2, f3, c3, f4)
#                 i += 1
#         print(i)
#         sys.exit()
            plotSeqDuplication([self.cdrSeqs['cdr1'].tolist(),
                              self.cdrSeqs['cdr2'].tolist(),
                              self.cdrSeqs['cdr3'].tolist(),
                              VH],
                             self.outputDir + self.name + 
                     '_Vdomain__Vdomain_ication.png',
                             ['CDR1', 'CDR2', 'CDR3', 'V Domain'],
                             'Duplication of V Domain Sequences')
            plotSeqDiversity([self.cdrSeqs['cdr1'].tolist(),
                          self.cdrSeqs['cdr2'].tolist(),
                          self.cdrSeqs['cdr3'].tolist(),
                          VH],
                         self.outputDir + self.name + 
                 '_Vdomain_diversity.png',
                         ['CDR1', 'CDR2', 'CDR3', 'V Domain'],
                         'Diversity of V Domain Sequences')
        gc.collect()
        
        plotSeqDuplication([self.cdrSeqs['cdr1'].tolist(),
                          self.cdrSeqs['cdr2'].tolist(),
                          self.cdrSeqs['cdr3'].tolist()
                          ],
                         self.outputDir + self.name + 
                 '_cdr_duplication.png',
                         ['CDR1', 'CDR2', 'CDR3'],
                         'Duplication of CDR Sequences')        
        
        plotSeqDiversity([self.cdrSeqs['cdr1'].tolist(),
                          self.cdrSeqs['cdr2'].tolist(),
                          self.cdrSeqs['cdr3'].tolist()
                        ],
                         self.outputDir + self.name + 
                 '_cdr_diversity.png',
                         ['CDR1', 'CDR2', 'CDR3'],
                         'Diversity of CDR Sequences')
        gc.collect()


    def analyze5UTR(self):
        print("The diversity of the upstream of IGV genes is being analyzed ... ")
#         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
#         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        self.upstreamFile = self.outputDir + self.name  
        self.upstreamFile += "_5utr_%.0f_%.0f.fasta" % (self.upstream[0],
                                                        self.upstream[1])
        self.alignInfoFile = self.outputDir + self.name
        self.alignInfoFile += "_align_info.tab"
        # extract upstream sequences 
        if (not exists(self.upstreamFile)):
            self.extractUpstreamSeqs()
        else:
            print("\tUpstream sequences file was found! ... " + self.upstreamFile)
        self.upstreamFile = os.path.abspath(self.upstreamFile)
        self.alignInfoFile = os.path.abspath(self.alignInfoFile)        
        # plot the distribution of sequence length
        expectLength = self.upstream[1] - self.upstream[0] + 1
        outputFile = self.upstreamFile.replace('.fasta', '_dist.png')
        plotSeqLenDist(self.upstreamFile, self.name, outputFile, self.format)
        outputFile = self.upstreamFile.replace('.fasta', '_dist_class.png')
        plotSeqLenDistClasses(self.upstreamFile, self.name, outputFile,
                              self.format)
        if (expectLength != Inf):       
            outputFile = self.upstreamFile.replace('.fasta', '_dist_short.png')
            plotSeqLenDist(self.upstreamFile, self.name, outputFile, self.format,
                           expectLength - 1)        
            outputFile = self.upstreamFile.replace('.fasta', '_dist_short_class.png')
            plotSeqLenDistClasses(self.upstreamFile, self.name, outputFile,
                              self.format, expectLength - 1)
            # # analyze intact secretion signals
            self.analyzeSequences(self.name, [expectLength, expectLength],
                                  startCodon=False, type='5utr', clusterMotifs=True)
        
    

    
