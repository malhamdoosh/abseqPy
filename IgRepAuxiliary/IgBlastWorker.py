'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

from multiprocessing import Process
from IgRepertoire.igRepUtils import runIgblastn, runIgblastp
from pandas.core.frame import DataFrame
import sys
import numpy as np


ANNOTATION_FIELDS = ['queryid', 'vgene', 'vqstart', 'vstart', 'vmismatches', 'vgaps',
              'identity', 'alignlen', 'bitscore',
                    'dgene', 'dqstart', 'dqend', 'dstart', 'dmismatches', 'dgaps', 
                    'jgene', 'jqstart', 'jqend', 'jstart', 'jmismatches', 'jgaps',  
                    'strand', 'stopcodon', 'v-jframe',
                    'fr1.start', 'fr1.end', 'fr1.mismatches', 'fr1.gaps',
                    'cdr1.start', 'cdr1.end', 'cdr1.mismatches', 'cdr1.gaps',
                    'fr2.start', 'fr2.end', 'fr2.mismatches', 'fr2.gaps',
                    'cdr2.start', 'cdr2.end', 'cdr2.mismatches', 'cdr2.gaps',
                    'fr3.start', 'fr3.end', 'fr3.mismatches', 'fr3.gaps',
                    'cdr3.start', 'cdr3.end', 'cdr3.mismatches', 'cdr3.gaps',
                    'fr4.start', 'fr4.end', 'fr4.mismatches', 'fr4.gaps'
                    ]

def getAnnotationFields(chain):
    if (chain == 'hv'):
        return ANNOTATION_FIELDS
    elif (chain in ['kv', 'lv']):
        return filter(lambda x: not x.startswith("d"), ANNOTATION_FIELDS)
    else:
        print('ERROR: unsupported chain type.')     
        sys.exit() 

def createCloneRecord(chain):
    cdrRecord = {}    
    for field in getAnnotationFields(chain):
        cdrRecord[field] =  np.nan
    return cdrRecord


def convertCloneRecordToOrderedList(cdrRecord, chain):    
    orderedList = []
    for field in getAnnotationFields(chain):
        orderedList.append(cdrRecord[field])
    
    return orderedList

def to_int(x):
    try:
        return int(x.strip())
    except:
        return None

def extractCDRInfo(blastOutput, chain):
    # Extract the top hits  
    print('\tExtracting top hit tables ... ' + blastOutput.split("/")[-1])
    # process igblast output and extract top hit 
    cloneAnnot = []
    filteredIDs = []
    line = ""
    
    warning = False
    with open(blastOutput) as blast:                   
        while(True):
            try:
                if (not line.startswith('# Query')): 
                    line = blast.readline()                       
                    if (not line):
                        break 
                    continue            
                cloneRecord = createCloneRecord(chain)
                cloneRecord['queryid'] = line.split()[2].strip()
                # parse  V-(D)-J rearrangement   
                line = blast.readline()
                while(line and 
                      not line.startswith('# Query') and
                       not line.startswith('# V-(D)-J rearrangement')):
                    line = blast.readline()
                if (not line):
                    filteredIDs.append(cloneRecord['queryid'])
                    break
                if (line.startswith('# Query')):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue
                line = blast.readline().strip().split('\t')
                cloneRecord['strand'] = 'forward' if line[-1] == '+' else 'reversed'
#                 print line, cloneRecord['strand']
#                 sys.exit()
                if (chain == 'hv'):
                    cloneRecord['stopcodon'] = line[4]
                    cloneRecord['v-jframe'] = line[5]
                    cloneRecord['vgene'] = line[0].split(',')[0]
                    cloneRecord['dgene'] = line[1].split(',')[0]
                    cloneRecord['jgene'] = line[2].split(',')[0]
                else:
                    cloneRecord['stopcodon'] = line[3]
                    cloneRecord['v-jframe'] = line[4]
                    cloneRecord['vgene'] = line[0].split(',')[0]                    
                    cloneRecord['jgene'] = line[1].split(',')[0]
                # parse Alignment Summary between query and top germline V gene
                line = ' '.join(line)
                while (line and 
                       not line.startswith('# Query') and
                       not line.startswith("# Alignment")):
                    line = blast.readline()
                if (not line):
                    filteredIDs.append(cloneRecord['queryid'])
                    break
                if (line.startswith('# Query')):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue
                line = blast.readline()
                for i in range(1,4):                
                    if (line.lower().startswith('fr' + `i`)):
                        line = line.split()
                        cloneRecord['fr%d.start' % i] = to_int(line[1])
                        cloneRecord['fr%d.end' % i] = to_int(line[2]) 
                        cloneRecord['fr%d.mismatches' % i] = to_int(line[5])
                        cloneRecord['fr%d.gaps' % i] = to_int(line[6])
                        line = blast.readline()
                    if (line.lower().startswith('cdr' + `i`)):
                        line = line.replace('(germline)', '').split()
                        cloneRecord['cdr%d.start' % i] = to_int(line[1])
                        cloneRecord['cdr%d.end' % i] = to_int(line[2]) 
                        cloneRecord['cdr%d.mismatches' % i] = to_int(line[5])
                        cloneRecord['cdr%d.gaps' % i] = to_int(line[6])
                        line = blast.readline()
                # parse alignment information between query and V, D and J genes
                while (line and 
                       not line.startswith('# Query') and
                       not line.startswith("# Fields")):
                    line = blast.readline()
                if (not line):
                    filteredIDs.append(cloneRecord['queryid'])
                    break
                if (line.startswith('# Query')):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue            
                line = blast.readline()
                noHits = to_int(line.split()[1])    
                if (noHits == 0):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue 
                # retrieve the top hit
                # parse the top V gene info
                line = blast.readline()
                if (not line.startswith("V")): 
                    filteredIDs.append(cloneRecord['queryid'])
                    continue
                hit = line.split()
                score = float(hit[-1]) 
                align = to_int(hit[4])   
                sStart = to_int(hit[10])                
                cloneRecord['identity'] = float(hit[3])
                cloneRecord['alignlen'] = align              
                cloneRecord['bitscore'] = score
                cloneRecord['vqstart'] = to_int(hit[8])
                cloneRecord['vstart'] = sStart
                cloneRecord['vmismatches'] = to_int(hit[5])
                cloneRecord['vgaps'] = to_int(hit[7])
                # parse the top D gene info
                line = blast.readline()
                while (line and
                       not line.startswith("# Query") and
                       not line.startswith("D") and 
                       not line.startswith("J")):
                    line = blast.readline()
                if (not line):
                    cloneAnnot.append(convertCloneRecordToOrderedList(cloneRecord, chain))
                    break
                if (line.startswith('# Query')):
                    cloneAnnot.append(convertCloneRecordToOrderedList(cloneRecord, chain))
                    continue         
                if (line.startswith("D")):
                    hit = line.split()
                    cloneRecord['dqstart'] = to_int(hit[8])
                    cloneRecord['dqend'] = to_int(hit[9])
                    cloneRecord['dstart'] = to_int(hit[10])
                    cloneRecord['dmismatches'] = to_int(hit[5])
                    cloneRecord['dgaps'] = to_int(hit[7])
                # parse the top J gene info
                while (line and
                   not line.startswith("# Query") and
                   not line.startswith("J")):
                    line = blast.readline()
                if (not line):
                    cloneAnnot.append(convertCloneRecordToOrderedList(cloneRecord, chain))
                    break
                if (line.startswith('# Query')):
                    cloneAnnot.append(convertCloneRecordToOrderedList(cloneRecord, chain))
                    continue 
                if (line.startswith("J")):
                    hit = line.split()
                    cloneRecord['jqstart'] = to_int(hit[8])
                    cloneRecord['jqend'] = to_int(hit[9])
                    cloneRecord['jstart'] = to_int(hit[10])
                    cloneRecord['jmismatches'] = to_int(hit[5])
                    cloneRecord['jgaps'] = to_int(hit[7])           
                cloneAnnot.append(convertCloneRecordToOrderedList(cloneRecord, chain)) 
            except Exception as e:                
#                 print(line, cloneRecord)
#                 raise e
                warning = True
                continue            
    if (len(cloneAnnot) > 0):
        # productive = no stop and in-frame
        # v-jframe: in-frame, out-of-frame, N/A (no J gene) 
        # stopcodon: yes, no
        cloneAnnot = DataFrame(cloneAnnot, columns =  getAnnotationFields(chain)) 
        cloneAnnot.index = cloneAnnot.queryid
        del cloneAnnot['queryid']
    else:
        cloneAnnot = DataFrame()
    if (warning):
        print("Warning: something went wrong while parsing %s" % (blastOutput))                                                            
    return (cloneAnnot, filteredIDs)


def analyzeSmallFile(fastaFile, chain, igBlastDB,
                     seqType='dna', threads=8): # , bitScore = 0
    # Run igblast
    if seqType.lower() == 'dna':
        blastOutput = runIgblastn(fastaFile, chain, threads, igBlastDB)   
    else:
        blastOutput = runIgblastp(fastaFile, chain, threads, igBlastDB)
    return extractCDRInfo(blastOutput, chain)
    

class IgBlastWorker(Process):
    def __init__(self, chain, igBlastDB, 
                 seqType, threads):
        super(IgBlastWorker, self).__init__() 
        self.chain = chain     
        self.igBlastDB = igBlastDB        
        self.seqType = seqType
        self.threads = threads
        self.tasksQueue = None
        self.resultsQueue = None
        self.exitQueue = None
#         self.args = None
    
#     def __init__(self, name, tasksQueue, resultsQueue):
#         self.name = name
#         self.tasksQueue = tasksQueue
#         self.resultsQueue = resultsQueue
    
    def run(self):
        while True:            
            nextTask = self.tasksQueue.get()
#             print("process has started a run... " + self.name)
            # poison pill check            
            if (nextTask is None):
                print("process has stopped ... " + self.name)
                self.exitQueue.put("exit")
#                 self.terminate()
                break
            try:
                result = analyzeSmallFile(nextTask, self.chain, self.igBlastDB,                                                                                                      
                                                     self.seqType, self.threads)                        
#                 print("process has completed analysis... " + self.name) 
                self.resultsQueue.put(result)            
            except Exception as e:
                print("An error occurred while processing " + nextTask.split('/')[-1])
                print(e)
                self.resultsQueue.put(None)
#                 raise
#                 sys.exit()
                continue                       
#             print("process has completed a run... " + self.name) 
        return
         