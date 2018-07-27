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

from abseqPy.IgMultiRepertoire.IgMultiRepertoire import IgMultiRepertoire
from abseqPy.argsParser import parseArgs
from abseqPy.config import VERSION
from abseqPy.utilities import PriorityPath
from abseqPy.logger import formattedTitle

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
                igRepertoire.start()

            print("The analysis started at " + startTimeStr)
            print("The analysis took {}".format(timedelta(seconds=int(round(time.time() - startTime)))))
            print("AbSeqPy version " + VERSION)
    except Exception as e:
        print("Unexpected error: " + str(e))
        print('-' * 60)
        traceback.print_exc(file=sys.stdout)
        print('-' * 60)
    finally:
        pass
# # Clean igblast files
#     os.system("rm " + blastOutput)
#     # Clean the output folder: remove fasta 
#     if (format == 'fastq'):
#         os.system("rm " + readFasta1)
#         os.system("rm " + readFasta2)
