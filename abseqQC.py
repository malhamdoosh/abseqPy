#!/usr/bin/env python

'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import sys
import os
from IgMultiRepertoire.IgMultiRepertoire import IgMultiRepertoire
import time
from argsParser import parseArgs
import traceback
from IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses
from config import VERSION
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.simplefilter(action = "ignore", category = DeprecationWarning)

__version__ = VERSION

# def fxn():
#     warnings.warn("deprecated", DeprecationWarning)
# 
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()

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
    logFile = None
    try:        
        argsVals = parseArgs()
        # print("Abseq output has been logged into " + argsVals.log)
        # logFile = open(argsVals.log, 'a')
        origStdout = sys.stdout
        # sys.stdout = logFile
        if (argsVals.task == 'all'):
            printFormattedTitle("Running the complete QC pipeline")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.runFastqc()
            igRepertoire.annotateClones(all=True)
            igRepertoire.analyzeAbundance(all=True)
            igRepertoire.analyzeProductivity(all=True)
            igRepertoire.analyzeDiversity(all=True)
        if (argsVals.task == 'fastqc'):
            printFormattedTitle("Sequencing QC Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.runFastqc()
        elif (argsVals.task == 'annotate'):
            printFormattedTitle("Clone Identification and Classification")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.annotateClones()
        elif (argsVals.task == 'abundance'):
            printFormattedTitle("IGV Abundance and QC Plots")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeAbundance() # estimateIGVDist()
        elif (argsVals.task == 'productivity'):
            printFormattedTitle("Clone Productivity Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeProductivity()
        elif (argsVals.task == 'diversity'):
            printFormattedTitle("Diversity Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeDiversity()
        elif (argsVals.task == 'secretion'):
            #analyze the sequences upstream of the IGV genes
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeSecretionSignal()
        elif (argsVals.task == '5utr'):
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyze5UTR()
        elif (argsVals.task == 'rsasimple'):
            printFormattedTitle("Simple Restriction Sites Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeRestrictionSitesSimple()
        elif (argsVals.task == 'rsa'):
            printFormattedTitle("Comprehensive Restriction Sites Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzeRestrictionSites()
        elif (argsVals.task == 'primer'):
            printFormattedTitle("Primer Specificity Analysis")
            igRepertoire = IgMultiRepertoire(argsVals)
            igRepertoire.analyzePrimerSpecificity()
        elif (argsVals.task == 'seqlen'):
            printFormattedTitle("Sequence Length Distribution")
            # calculate the distribution of sequence lengths of a sample
            argsVals.outdir += 'annot/'
            os.system("mkdir " + argsVals.outdir)
            outputFile =  argsVals.outdir + argsVals.name + '_seq_length_dist.png'
            plotSeqLenDist(argsVals.f1, argsVals.name, outputFile, argsVals.fmt,
                           maxbins=-1)
        elif (argsVals.task == 'seqlenclass'):
            # calculate the distribution of sequences in different IGV families
            # input file must be a file of IGV genes
            argsVals.o += 'annot/'
            os.system("mkdir " + argsVals.outdir)
            outputFile =  argsVals.outdir + argsVals.name + '_length_dist_classes.png'
            plotSeqLenDistClasses(argsVals.f1, argsVals.name, outputFile, argsVals.fmt)

        print ("The analysis started at " + startTimeStr)
        print "The analysis took %.2f  minutes!!" % ((time.time() - t) / 60)
        print("Abseq Version " + VERSION)

    except Exception as e:
        print("Unexpected error: " + str(e))
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
    finally:        
        if logFile is not None:
            logFile.close()
        
    
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


