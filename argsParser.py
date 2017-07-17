'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 
import os
from numpy import Inf
from os.path import abspath

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

PROGRAM_VALID_ARGS = ['-task', '-chain', '-name',
                 '-f1', '-f2', '-fmt', '-o', '-merge', '-merger',  
                 '-seqtype', '-threads', '-db',  
                  '-bitscore', '-alignlen', '-sstart', '-actualqstart',     
                  '-trim5' ,'-trim3',  '-fr4cut',            
                   '-sites',
                  '-primer', 
                  '-5end', '-3end', 
                  '-5endoffset',
                 '-upstream',
                 '-report-interim'
                  ]

def parseArgs(args):    
    argVals = {}
    argVals["cmd"] = ' '.join(args)
    for i in range(1, len(args), 2):
#         print(args[i].lower(), PROGRAM_VALID_ARGS)
        if (args[i].lower() in PROGRAM_VALID_ARGS):
#             print(args[i], args[i+1])
            argVals[args[i][1:].lower()] = args[i+1]
        else:
            print(args[i], PROGRAM_VALID_ARGS)
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
        
    # input sequencing files
    if (argVals.get('f1', None) is None):
        raise Exception("One sequence file at least must be provided.")
    elif not os.path.exists(argVals['f1']):        
        raise Exception("-f1 file not found!")    
    else:
        argVals['f1'] = abspath(argVals['f1'])    
        
    if (argVals.get('f2', None) is None):
        argVals['f2'] = None
    elif not os.path.exists(argVals['f2']):
        raise Exception("-f2 file not found!")
    else:
        argVals['f2'] = abspath(argVals['f2'])
        
    # merging options  
    if (argVals.get('merge', None) is None):
        argVals['merge'] = 'no'
    else:
        argVals['merge'] = argVals['merge'].lower()        
    if (argVals['merge'] == 'yes' and (argVals['f2'] is None)):
        raise Exception("The merger requires two sequence files (use both -f1 and -f2).")
    else:
        if (argVals.get('merger', None) is None):
            argVals['merger'] = 'flash'
        else:
            argVals['merger'] = argVals['merger'].lower()
  
    ## output directory and parse sample name 
    if (argVals.get('o', None) is None):
        argVals['o'] = './'     
    if (argVals.get('name', None) is not None):  
        argVals['o'] += "/" + argVals['name']
    else:
        f1name = argVals['f1'].split('/')[-1]
        if (f1name.find("_R") != -1 and 
            (argVals['merge'] == 'yes' or argVals['task'] == "fastqc")):
            ext = '_' + f1name.split('_')[-1]
        else:
            ext = f1name[f1name.find('.'):]    
        argVals['o'] += "/"+f1name.replace(ext, '')
        sampleName = f1name.split("_")[0] + '_'  
        sampleName += f1name.split("_")[-1].split('.')[0]
        argVals['name'] = sampleName
    argVals['o'] = (abspath(argVals['o']) + '/').replace('//', '/')
    os.system("mkdir -p " + argVals['o'])
    
        
    
    # (argVals['merge'] == 'no' and noFiles == 2) or 
    
        
    if (argVals.get('fmt', None) is None):
        argVals['fmt'] = 'fastq'
    else:
        argVals['fmt'] = argVals['fmt'].lower()
    
    if (argVals.get('bitscore', None) is None):
        argVals['bitscore'] = [0, Inf]
    else:
        argVals['bitscore'] = extractRanges(argVals['bitscore'])[0]
        
    if (argVals['task'] in ['secretion', '5utr']):
        if (argVals.get('upstream', None) is None):
            argVals['upstream'] = [1, Inf] 
        else:
            argVals['upstream'] = extractRanges(argVals['upstream'], 1)[0]
    if (argVals['task'] in ['rsa', 'rsasimple']):
        if (argVals.get('sites', None) is None):
            raise Exception("Restriction sites should be provided.")            
        argVals['sites'] = abspath(argVals['sites'])
    if (argVals['task'] in ['diversity', 'productivity', 'all']):  
        # actualqstart is a 1-based index      
        if argVals.get('actualqstart', None) is not None:
            argVals['actualqstart'] = int(argVals['actualqstart']) - 1
        else:
            argVals['actualqstart'] = -1        
        if (argVals.get('fr4cut', None) is None):
            argVals['fr4cut'] = False
        else:
            argVals['fr4cut'] = (argVals['fr4cut'].upper() == 'YES')    
    if (argVals.get('trim5', None) is not None):
        argVals['trim5'] = int(argVals['trim5']) - 1
    else:
        argVals['trim5'] = 0
    if (argVals.get('trim3', None) is not None):
        argVals['trim3'] = int(argVals['trim3'])
    else:
        argVals['trim3'] = 0
    if (argVals['task'] == 'primer'):
        argVals['5end'] = argVals.get('5end', None)
        if argVals['5end']:
            argVals['5end'] = abspath(argVals['5end'])
        argVals['3end'] = argVals.get('3end', None)
        if argVals['3end']:
            argVals['3end'] = abspath(argVals['3end'])
        if (argVals.get('5endoffset', None) is None):
            argVals['5endoffset'] = 0
        else:
            argVals['5endoffset'] = int(argVals['5endoffset'])            
        
    if (argVals.get('sstart', None) is None):
        argVals['sstart'] = [1, Inf]
    else:
        argVals['sstart'] = extractRanges(argVals['sstart'])[0]  
    
    if (argVals.get('alignlen', None) is None):
        argVals['alignlen'] = [0, Inf]
    else:
        argVals['alignlen'] = extractRanges(argVals['alignlen'])[0]    
            
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
        
    if (argVals.get('report-interim', None) is None):
        argVals['report-interim'] = False
    else:
        argVals['report-interim'] = (argVals['report-interim'].upper() == 'YES')    
    argVals['log'] = argVals['o'] + argVals['name'] + ".log"    
        
    return argVals





