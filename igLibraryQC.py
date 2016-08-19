#!/usr/bin/env python

## Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)

import sys
import os
from numpy import Inf
from IgRepertoire.IgRepertoire import IgRepertoire
import time
from os.path import abspath
import traceback
from IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses
from argsParser import parseArgs
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

def printFormattedTitle(title):
    print "-" * 100
    print "|" + " " * 98 + "|"
    print "|" + " " * ((98 - len(title)) / 2) + title + " " * ((98 - len(title))/2 + (98 - len(title))%2) + "|"
    print "|" + " " * 98 + "|"
    print "-" * 100
    sys.stdout.flush()

def main():
    
    startTimeStr = time.strftime("%Y-%m-%d %H:%M:%S")
    t = time.time()
    try:
        argsVals = parseArgs(sys.argv)
        #print(argsVals)        
        if (argsVals['task'] == 'fastqc'): 
            printFormattedTitle("Sequencing QC Analysis")
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.runFastqc()
        elif (argsVals['task'] == 'annotate'): 
            printFormattedTitle("Clone Identification and Classification")
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.annotateClones()
        elif (argsVals['task'] == 'abundance'): 
            printFormattedTitle("IGV Abundance and QC Plots")
            igRepertoire = IgRepertoire(argsVals)
            igRepertoire.analyzeAbundance() # estimateIGVDist()  
        elif (argsVals['task'] == 'productivity'):
            printFormattedTitle("Clone Productivity Analysis")
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.analyzeProductivity()  
        elif (argsVals['task'] == 'diversity'):
            printFormattedTitle("CDR Sequence Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeDiversity()       
        elif (argsVals['task'] == 'secretion'):
            #analyze the sequences upstream of the IGV genes
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.analyzeSecretionSignal()
        elif (argsVals['task'] == '5utr'):
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.analyze5UTR()        
        elif (argsVals['task'] == 'enzymesimple'):
            printFormattedTitle("Simple Restriction Sites Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeRestrictionSitesSimple()
        elif (argsVals['task'] == 'enzymes'):
            printFormattedTitle("Comprehensive Restriction Sites Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeRestrictionSites()
        elif (argsVals['task'] == 'primer'):
            printFormattedTitle("Primer Specificity Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzePrimerSpecificity()
        elif (argsVals['task'] == 'seqlen'):
            printFormattedTitle("Sequence Length Distribution")
            # calculate the distribution of sequence lengths of a sample
            argsVals['o'] += 'annot/'
            os.system("mkdir " + argsVals['o'])
            outputFile =  argsVals['o'] + argsVals['name'] + '_seq_length_dist.png'            
            plotSeqLenDist(argsVals['f1'], argsVals['name'], outputFile, argsVals['fmt'],
                           maxbins=-1)
        elif (argsVals['task'] == 'seqlenclass'):
            # calculate the distribution of sequences in different IGV families
            # input file must be a file of IGV genes
            argsVals['o'] += 'annot/'
            os.system("mkdir " + argsVals['o'])
            outputFile =  argsVals['o'] + argsVals['name'] + '_length_dist_classes.png'
            plotSeqLenDistClasses(argsVals['f1'], argsVals['name'], outputFile, argsVals['fmt'])
        print ("The analysis started at " + startTimeStr)
        print "The analysis took %.2f  minutes!!" % ((time.time() - t) / 60)
    
    except Exception as e:
        print("Unexpected error: " + str(e))
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
        
    
#TODO: generate HTML report for each analysis and give name to the report
'''
Report parameters
Organize figures 
'''
            
# # Clean igblast files
#     os.system("rm " + blastOutput)
#     # Clean the output folder: remove fasta 
#     if (format == 'fastq'):
#         os.system("rm " + readFasta1)
#         os.system("rm " + readFasta2)


