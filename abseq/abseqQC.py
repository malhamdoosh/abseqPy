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

from datetime import timedelta

from abseq.IgMultiRepertoire.IgMultiRepertoire import IgMultiRepertoire
from abseq.argsParser import parseArgs
from abseq.config import VERSION
from abseq.utilities import PriorityPath
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
    startTime = time.time()
    try:
        with PriorityPath():
            argsVals = parseArgs()

            with IgMultiRepertoire(argsVals) as igRepertoire:
                # show a pretty banner before beginning analysis
                print(formattedTitle(argsVals.task, argsVals.yaml is not None))
                igRepertoire.rockNRoll()

            print("The analysis started at " + startTimeStr)
            print("The analysis took {}".format(timedelta(seconds=int(round(time.time() - startTime)))))
            print("Abseq Version " + VERSION)
    except Exception as e:
        print("Unexpected error: " + str(e))
        print('-' * 60)
        traceback.print_exc(file=sys.stdout)
        print('-' * 60)
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
