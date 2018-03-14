#!/usr/bin/env python

'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import sys
import time
import traceback
import warnings

from abseq.IgMultiRepertoire.IgMultiRepertoire import IgMultiRepertoire
from abseq.argsParser import parseArgs
from abseq.config import VERSION, PriorityPath
from abseq.logger import formattedTitle

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)

__version__ = VERSION


# def fxn():
#     warnings.warn("deprecated", DeprecationWarning)
# 
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()


def main():
    startTimeStr = time.strftime("%Y-%m-%d %H:%M:%S")
    t = time.time()
    try:
        with PriorityPath():
            argsVals = parseArgs()
            with IgMultiRepertoire(argsVals) as igRepertoire:

                # show a pretty banner before beginning analysis
                print(formattedTitle(argsVals.task))

                if argsVals.task == 'all':
                    igRepertoire.runFastqc()
                    igRepertoire.annotateClones(all=True)
                    igRepertoire.analyzeAbundance(all=True)
                    igRepertoire.analyzeProductivity(all=True)
                    igRepertoire.analyzeDiversity(all=True)
                if argsVals.task == 'fastqc':
                    igRepertoire.runFastqc()
                elif argsVals.task == 'annotate':
                    igRepertoire.annotateClones()
                elif argsVals.task == 'abundance':
                    igRepertoire.analyzeAbundance()  # estimateIGVDist()
                elif argsVals.task == 'productivity':
                    igRepertoire.analyzeProductivity()
                elif argsVals.task == 'diversity':
                    igRepertoire.analyzeDiversity()
                elif argsVals.task == 'secretion':
                    # analyze the sequences upstream of the IGV genes
                    igRepertoire.analyzeSecretionSignal()
                elif argsVals.task == '5utr':
                    igRepertoire.analyze5UTR()
                elif argsVals.task == 'rsasimple':
                    igRepertoire.analyzeRestrictionSitesSimple()
                elif argsVals.task == 'rsa':
                    igRepertoire.analyzeRestrictionSites()
                elif argsVals.task == 'primer':
                    igRepertoire.analyzePrimerSpecificity()
                elif argsVals.task == 'seqlen':
                    # calculate the distribution of sequence lengths of a sample
                    igRepertoire.analyzeSeqLen()
                elif argsVals.task == 'seqlenclass':
                    # calculate the distribution of sequences in different IGV families
                    # input file must be a file of IGV genes
                    igRepertoire.analyzeSeqLen(klass=True)

            print("The analysis started at " + startTimeStr)
            print("The analysis took %.2f  minutes!!" % ((time.time() - t) / 60))
            print("Abseq Version " + VERSION)
    except Exception as e:
        print("Unexpected error: " + str(e))
        print '-' * 60
        traceback.print_exc(file=sys.stdout)
        print '-' * 60
    finally:
        pass

# TODO: generate HTML report for each analysis and give name to the report
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
