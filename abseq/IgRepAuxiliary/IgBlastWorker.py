'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import os
import numpy as np

from multiprocessing import Process
from pandas.core.frame import DataFrame

from abseq.IgRepertoire.igRepUtils import runIgblastn, runIgblastp
from abseq.logger import printto, LEVEL

ANNOTATION_FIELDS = ['queryid', 'vgene', 'vqstart', 'vstart', 'vmismatches', 'vgaps',
                    'identity', 'alignlen', 'bitscore',
                    'dgene', 'dqstart', 'dqend', 'dstart', 'dmismatches', 'dgaps', 
                    'jgene', 'jqstart', 'jqend', 'jstart', 'jmismatches', 'jgaps',  
                    'strand', 'stopcodon', 'v-jframe',
                    'fr1.start', 'fr1.end', 'fr1.mismatches', 'fr1.gaps',
                    'cdr1.start', 'cdr1.end', 'cdr1.mismatches', 'cdr1.gaps',
                    'fr2.start', 'fr2.end', 'fr2.mismatches', 'fr2.gaps',
                    'cdr2.start', 'cdr2.end', 'cdr2.mismatches', 'cdr2.gaps',
                    'fr3.start', 'fr3g.end', 'fr3g.mismatches', 'fr3g.gaps', 'fr3.end',
                    'cdr3g.start', 'cdr3g.end', 'cdr3g.mismatches', 'cdr3g.gaps', 'cdr3.start', 'cdr3.end',
                    'fr4.start', 'fr4.end', 'fr4.mismatches', 'fr4.gaps'
                    ]


def getAnnotationFields(chain):
    if chain == 'hv':
        return ANNOTATION_FIELDS
    elif chain in ['kv', 'lv']:
        return filter(lambda x: not x.startswith("d"), ANNOTATION_FIELDS)
    else:
        # should never happen (argparse takes care of this for us)
        raise Exception("Unsupported chain type")


def createCloneRecord(chain):
    cdrRecord = {}    
    for field in getAnnotationFields(chain):
        cdrRecord[field] = np.nan
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


def extractCDRInfo(blastOutput, chain, stream=None):
    # Extract the top hits  
    printto(stream, '\tExtracting top hit tables ... ' + os.path.basename(blastOutput))
    # process igblast output and extract top hit 
    cloneAnnot = []
    filteredIDs = []
    line = ""
    
    warning = False
    # RE: parsing IGBLAST:
    # VDJ junction details MAY give N/A instead of just missing:
    # eg:
    # V-(D)-J junction details based on top germline gene matches
    # (V end, V-D junction, D region, D-J junction, J start).  Note that possible
    # overlapping nucleotides at VDJ junction (i.e, nucleotides that could
    # be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are
    # not included under the V, D, or J gene itself
    # GGGTC  TGTTCACGAGGGCATCTGTGTCCTGTTTTTAGGTTCTCCTCCC  TTTTGAC  N/A  N/A

    # it also has a variable number of hits (depending on presence of region)
    # EG:
    # V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
    # CCTCT  N/A  GGTGT

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

                # XXX: the or len(line) == 8 may happen to light chains too, when there is
                # a rogue D-gene that was a hit. It then follows heavy chain's indexing
                if (chain == 'hv') or len(line) == 8:
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
                line = ' '.join(line)

                # Parse Sub-region analysis and ignore it if there's no CDR3 hit by IGBLAST
                while line and \
                        not line.startswith("# Alignment") and \
                        not line.startswith("# Sub-region") and \
                        not line.startswith("# Query"):
                    line = blast.readline()

                # EOF
                if not line:
                    filteredIDs.append(cloneRecord['queryid'])
                    break

                # there's no # Sub-region, nor is there # Alignment.
                if line.startswith("# Query"):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue

                # this implies that IGBLAST successfully classified a CDR3 sequence
                if line.startswith("# Sub-region"):
                    line = blast.readline()
                    subregionData = line.split()
                    assert subregionData[0] == 'CDR3'
                    if len(subregionData) >= 3 and subregionData[-1].isdigit() and subregionData[-2].isdigit():
                        cloneRecord['cdr3.start'] = to_int(subregionData[-2])
                        cloneRecord['cdr3.end'] = to_int(subregionData[-1])
                        # true FR3 end is at position cdr3.start - 1
                        # (the alignment table only tells us the FR3 germline)
                        # but since fr3.start always begins in the germline, there's no special field for that
                        # and is assumed that fr3.start == fr3g.start
                        cloneRecord['fr3.end'] = cloneRecord['cdr3.start'] - 1

                # parse Alignment Summary between query and top germline V gene
                while (line and
                       not line.startswith('# Query') and
                       not line.startswith("# Alignment")):
                    line = blast.readline()
                if not line:
                    filteredIDs.append(cloneRecord['queryid'])
                    break
                if line.startswith('# Query'):
                    filteredIDs.append(cloneRecord['queryid'])
                    continue
                line = blast.readline()
                for i in range(1, 4):
                    if line.lower().startswith('fr' + str(i)):
                        line = line.split()
                        cloneRecord['fr%d.start' % i] = to_int(line[1])
                        cloneRecord['fr%d%s.end' % (i, 'g' if i == 3 else '')] = to_int(line[2])
                        cloneRecord['fr%d%s.mismatches' % (i, 'g' if i == 3 else '')] = to_int(line[5])
                        cloneRecord['fr%d%s.gaps' % (i, 'g' if i == 3 else '')] = to_int(line[6])
                        line = blast.readline()
                    if line.lower().startswith('cdr' + str(i)):
                        line = line.replace('(germline)', '').split()
                        cloneRecord['cdr%d%s.start' % (i, 'g' if i == 3 else '')] = to_int(line[1])
                        cloneRecord['cdr%d%s.end' % (i, 'g' if i == 3 else '')] = to_int(line[2])
                        cloneRecord['cdr%d%s.mismatches' % (i, 'g' if i == 3 else '')] = to_int(line[5])
                        cloneRecord['cdr%d%s.gaps' % (i, 'g' if i == 3 else '')] = to_int(line[6])
                        line = blast.readline()

                # if the CDR3 region wasn't identified by IgBlast, we can't get FR3 end, so we fallback to
                # FR3 germline end. Since CDR3.start and CDR3.end isn't really used unitl the refinement process,
                # we ignore the fallback options for them.
                if np.isnan(cloneRecord['fr3.end']):
                    cloneRecord['fr3.end'] = cloneRecord['fr3g.end']

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
        printto(stream, "Warning: something went wrong while parsing %s" % (blastOutput), LEVEL.WARN)
    return (cloneAnnot, filteredIDs)


def analyzeSmallFile(fastaFile, chain, igBlastDB,
                     seqType='dna', threads=8, outdir="", stream=None): # , bitScore = 0
    # Run igblast
    if seqType.lower() == 'dna':
        blastOutput = runIgblastn(fastaFile, chain, threads, igBlastDB, outputDir=outdir, stream=stream)
    else:
        # argparse already checks that it's ether dna or protein, so nothing fishy can pass into else statement here
        blastOutput = runIgblastp(fastaFile, chain, threads, igBlastDB, outputDir=outdir, stream=stream)
    return extractCDRInfo(blastOutput, chain, stream=stream)
    

class IgBlastWorker(Process):
    def __init__(self, chain, igBlastDB, 
                 seqType, threads, stream=None):
        super(IgBlastWorker, self).__init__() 
        self.chain = chain     
        self.igBlastDB = igBlastDB        
        self.seqType = seqType
        self.threads = threads
        self.tasksQueue = None
        self.resultsQueue = None
        self.exitQueue = None
        self.stream = stream
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
            if nextTask is None:
                printto(self.stream, "process has stopped ... " + self.name)
                self.exitQueue.put("exit")
#                 self.terminate()
                break
            try:
                result = analyzeSmallFile(nextTask, self.chain, self.igBlastDB,                                                                                                      
                                          self.seqType, self.threads, stream=self.stream)
#                 print("process has completed analysis... " + self.name) 
                self.resultsQueue.put(result)            
            except Exception as e:
                printto(self.stream, "An error occurred while processing " + os.path.basename(nextTask), LEVEL.EXCEPT)
                self.resultsQueue.put(None)
#                 raise
#                 sys.exit()
                continue                       
#             print("process has completed a run... " + self.name) 
        return
