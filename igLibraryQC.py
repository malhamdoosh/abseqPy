#!/usr/bin/env python

## Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)

import sys
import os
from numpy import Inf
from IgRepertoire import IgRepertoire
from NGSutils import plotSeqLenDist
from NGSutils import plotSeqLenDistClasses
import time
from os.path import abspath
from sys import argv
import traceback
  

def extractRanges(strRanges, expNoRanges=2):
    numRanges = []
    ranges = strRanges.split(',')
    if (len(ranges) > expNoRanges):
        raise Exception("Number of bitScore, alignLen and sstart ranges should match the number of files")
        
    for i in range(len(ranges)):
        scores = ranges[i].split('-')
        if (len(scores) == 2):
            numRanges.append([float(scores[0]), float(scores[1])])
            if (numRanges[-1][0] >= numRanges[-1][1]):
                raise Exception("Invalid ranges " + strRanges)
                
        else:
            numRanges.append([float(scores[0]), Inf])
    if (len(numRanges) == 1 < expNoRanges):
        numRanges = numRanges * expNoRanges
    return numRanges

def parseArgs(args):
    validArgs = ['-task', '-chain', '-name',
                 '-f1', '-f2', '-o', '-merge', '-fmt', '-seqtype',
                  '-bitscore', '-primer', '-alignlen', '-sstart',
                  '-db', '-threads', '-merger', '-upstream',
                  '-sites', '-5end', '-3end', '-actualqstart',
                  '-5endoffset', '-fr4cut']
    argVals = {}
    for i in range(1, len(args), 2):
#         print(args[i].lower(), validArgs)
        if (args[i].lower() in validArgs):
#             print(args[i], args[i+1])
            argVals[args[i].replace('-', '').lower()] = args[i+1]
        else:
            raise Exception(args[i] + ' is invalid argument.')
    
    if (argVals.get('chain', None) is None):
        argVals['chain'] = "hv"
    else:
        argVals['chain'] = argVals['chain'].lower()
    if (argVals.get('task', None) is None):
        argVals['task'] = 'abundance'
    else:
        argVals['task'] = argVals['task'].lower()
        
    if (argVals.get('seqtype', None) is None):
        argVals['seqtype'] = 'dna'
        
    if (argVals.get('f1', None) is None):
        raise Exception("One sequence file at least must be provided.")    
    argVals['f1'] = abspath(argVals['f1'])
        
    if (argVals.get('f2', None) is not None):
        noFiles = 2
        argVals['f2'] = abspath(argVals['f2'])
    else:
        argVals['f2'] = None
        noFiles = 1
                
    if (argVals.get('o', None) is None):
        argVals['o'] = './'
    f1name = argVals['f1'].split('/')[-1]
    if (f1name.find("_R") != -1):
        ext = '_' + f1name.split('_')[-1]
    else:
        ext = f1name[f1name.find('.'):]    
    argVals['o'] += "/"+f1name.replace(ext, '')
    argVals['o'] = (abspath(argVals['o']) + '/').replace('//', '/')
    os.system("mkdir -p " + argVals['o'])
    
    if (argVals.get('merge', None) is None):
        argVals['merge'] = 'no'
    else:
        argVals['merge'] = argVals['merge'].lower()
        
    if (argVals['merge'] == 'yes' and (argVals['f2'] is None)):
        raise Exception("The merger requires two sequence files (f1 and f2).")
    
    if (argVals.get('merger', None) is None):
        argVals['merger'] = 'flash'
    else:
        argVals['merger'] = argVals['merger'].lower()
    
    if ((argVals['merge'] == 'no' and noFiles == 2) or argVals.get('name', None) is None):
        sampleName = argVals['f1'].split('/')[-1].split("_")[0] + '_'  
        sampleName += argVals['f1'].split('/')[-1].split("_")[-1].split('.')[0]
        argVals['name'] = sampleName
        
    if (argVals.get('fmt', None) is None):
        argVals['fmt'] = 'fastq'
    else:
        argVals['fmt'] = argVals['fmt'].lower()
    
    if (argVals.get('bitscore', None) is None):
        argVals['bitscore'] = [[0, Inf]] * noFiles
    else:
        argVals['bitscore'] = extractRanges(argVals['bitscore'], noFiles)
        
    if (argVals['task'] in ['upighv', '5utr']):
        if (argVals.get('upstream', None) is None):
            argVals['upstream'] = [1, Inf] 
        else:
            argVals['upstream'] = extractRanges(argVals['upstream'], 1)[0]
    if (argVals['task'] in ['enzymes', 'enzymesimple']):
        if (argVals.get('sites', None) is None):
            raise Exception("Restriction sites should be provided.")            
        argVals['sites'] = abspath(argVals['sites'])
    if (argVals['task'] == 'cdrs'):
        argVals['5end'] = argVals.get('5end', None)
        if argVals['5end']:
            argVals['5end'] = abspath(argVals['5end'])
        argVals['3end'] = argVals.get('3end', None)
        if argVals['3end']:
            argVals['3end'] = abspath(argVals['3end'])
        if argVals.get('actualqstart', None) is not None:
            argVals['actualqstart'] = int(argVals['actualqstart']) - 1
        else:
            argVals['actualqstart'] = -1
        if (argVals.get('5endoffset', None) is None):
            argVals['5endoffset'] = 0
        else:
            argVals['5endoffset'] = int(argVals['5endoffset'])
        if (argVals.get('fr4cut', None) is None):
            argVals['fr4cut'] = False
        else:
            argVals['fr4cut'] = (argVals['fr4cut'].upper() == 'YES')
        
    if (argVals.get('sstart', None) is None):
        argVals['sstart'] = [[1, Inf]] * noFiles
    else:
        argVals['sstart'] = extractRanges(argVals['sstart'], noFiles)        
    
    if (argVals.get('alignlen', None) is None):
        argVals['alignlen'] = [[0, Inf]] * noFiles
    else:
        argVals['alignlen'] = extractRanges(argVals['alignlen'], noFiles)        
            
    if (argVals.get('primer', None) is None):
        argVals['primer'] = -1
    else:
        argVals['primer'] = int(argVals['primer'])
        
    if (argVals.get('db', None) is None):
        argVals['db'] = '$IGBLASTDB'
    else:
        argVals['db'] = os.path.abspath(argVals['db'])
        
    if (argVals.get('threads', None) is None):
        argVals['threads'] = 8
    else:
        argVals['threads'] = int(argVals['threads'])
        
    return argVals

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
         
        if (argsVals['task'] == 'abundance'): 
            printFormattedTitle("IGV Abundance and QC Plots")
            igRepertoire = IgRepertoire(argsVals)
            igRepertoire.analyzeAbundance() # estimateIGVDist()
        elif (argsVals['task'] == 'seqlen'):
            printFormattedTitle("Sequence Length Distribution")
            # calculate the distribution of sequence lengths of a sample
            argsVals['o'] += 'abundance/'
            os.system("mkdir " + argsVals['o'])
            outputFile =  argsVals['o'] + argsVals['name'] + '_seq_length_dist.png'            
            plotSeqLenDist(argsVals['f1'], argsVals['name'], outputFile, argsVals['fmt'],
                           maxbins=-1)
        elif (argsVals['task'] == 'seqlenclass'):
            # calculate the distribution of sequences in different IGV families
            # input file must be a file of IGV genes
            argsVals['o'] += 'abundance/'
            os.system("mkdir " + argsVals['o'])
            outputFile =  argsVals['o'] + argsVals['name'] + '_length_dist_classes.png'
            plotSeqLenDistClasses(argsVals['f1'], argsVals['name'], outputFile, argsVals['fmt'])
        elif (argsVals['task'] == 'upighv'):
            #analyze the sequences upstream of the IGV genes
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.analyzeSecretionSignal()
        elif (argsVals['task'] == '5utr'):
            igRepertoire = IgRepertoire(argsVals)        
            igRepertoire.analyze5UTR()
        elif (argsVals['task'] == 'cdrs'):
            printFormattedTitle("CDR Sequence Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeCDRs()
        elif (argsVals['task'] == 'enzymesimple'):
            printFormattedTitle("Simple Restriction Sites Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeRestrictionSitesSimple()
        elif (argsVals['task'] == 'enzymes'):
            printFormattedTitle("Comprehensive Restriction Sites Analysis")
            igRepertoire = IgRepertoire(argsVals)    
            igRepertoire.analyzeRestrictionSites()
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


