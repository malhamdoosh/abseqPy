'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import re
import sys

from numpy import isnan, nan
from multiprocessing import Queue, Manager
from math import ceil
from pandas.core.frame import DataFrame

from abseqPy.IgRepAuxiliary.seqUtils import readSeqFileIntoDict
from abseqPy.utilities import hasLargeMem
from abseqPy.IgRepAuxiliary.RestrictionSitesScanner import RestrictionSitesScanner
from abseqPy.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseqPy.logger import printto, LEVEL


def initSimpleRSAStats(sites):
    # replace all with stats dict to be passed to a reporting function
    stats = {"siteHitsCount": {}, "siteHitSeqsCount": {}, "siteHitsSeqsIDs": {}, "seqsCutByAny": 0, "total": 0}
    for site in sites.keys():
        stats["siteHitsCount"][site] = 0
        stats["siteHitSeqsCount"][site] = 0
        stats["siteHitsSeqsIDs"][site] = []
    return stats
        

def postProcessRSA(stats, sitesInfo, stream=None):
    rsaResults = []
    sites = sorted(stats["siteHitSeqsCount"], key=stats["siteHitSeqsCount"].get)    
    for site in sites:
        rsaResults.append([site, sitesInfo[site],
                           stats["siteHitsCount"][site],
                           stats["siteHitsCount"][site] * 100.0 / sum(stats["siteHitsCount"].values()),
                           stats["siteHitSeqsCount"][site],
                           stats["siteHitSeqsCount"][site] * 100.0 / stats["total"]])
    rsaResults.append(["Cut by any", nan, nan, nan, stats["seqsCutByAny"],
                       stats["seqsCutByAny"] * 100.0 / stats["total"]])
    rsaResults.append(["Total", nan, nan, nan, stats["total"], 100])
    rsaResults = DataFrame(rsaResults, columns=["Enzyme", "Restriction Site", "No.Hits", "Percentage of Hits (%)",
                                                "No.Molecules", "Percentage of Molecules (%)"])
    overlapResults = {"order1": stats["siteHitsSeqsIDs"]}
    if len(stats["siteHitsSeqsIDs"]) >= 3:
        overlapResults["order2"] = calcRSAOverlapOrder2(stats["siteHitsSeqsIDs"], sites, stream=stream)
    return rsaResults, overlapResults


def calcRSAOverlapOrder2(order1, sites, stream=None):
    printto(stream, "The 2nd order overlapping matrix is being calculated using Jaccard Index ... ")
    overlap = []
    for site1 in sites:
        overlap.append([])
        for site2 in sites:
            inter = len(order1[site1].intersection(order1[site2])) * 1.0
            uni = len(order1[site1].union(order1[site2]))
            if uni != 0:
                overlap[-1].append(inter / uni)
            else:
                overlap[-1].append(1)
    overlap = DataFrame(overlap, columns=sites, index=sites)
    # overlap = linkage(overlap)
    return overlap


def scanRestrictionSitesSimple(name, readFile, format, cloneAnnot, sitesFile, threads, stream=None):
    sitesInfo = loadRestrictionSites(sitesFile, stream=stream)
    seqsPerWorker = len(sitesInfo)
    workers = []   
    try:
        m = Manager()        
        records = m.dict()
        readSeqFileIntoDict(readFile, format, records, stream=stream)
        queryIds = cloneAnnot.index
        noSeqs = len(queryIds)
        printto(stream, "{:,} restriction sites are being scanned for {:,} sequences ..."
                .format(len(sitesInfo), noSeqs))
        # Parallel implementation of the refinement
        totalTasks = int(ceil(noSeqs * 1.0 / seqsPerWorker)) 
        tasks = Queue()      
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs, desc="sequences", stream=stream)
        threads = min(threads, totalTasks)

        # Initialize workers
        for i in range(threads):
            w = RestrictionSitesScanner(records, cloneAnnot, procCounter,
                                        sitesInfo.copy(), simpleScan=True, stream=stream)
            w.tasksQueue = tasks
            w.exitQueue = exitQueue  
            w.resultsQueue = resultsQueue    
            workers.append(w)
            w.start()

        assert totalTasks > 0

        for i in range(totalTasks):
            ids = queryIds[i * seqsPerWorker:(i+1) * seqsPerWorker]
            tasks.put(ids)

        # Add a poison pill for each worker
        for i in range(threads + 10):
            tasks.put(None)

        # Wait all process workers to terminate
        i = 0
        while i < threads:    
            m = exitQueue.get()
            if m == "exit":
                i += 1

        # Collect results
        printto(stream, "All workers have completed their tasks successfully.")
        printto(stream, "Results are being collated from all workers ...")
        stats = collectRSASimpleResults(sitesInfo, resultsQueue, totalTasks, noSeqs, stream=stream)
        (rsaResults, overlapResults) = postProcessRSA(stats, sitesInfo, stream=stream)
        # End of parallel implementation
        printto(stream, "Results were collated successfully.")

    except Exception as e:
        printto(stream, "Something went wrong during the RSA scanning process!")
        raise e
    finally:        
        for w in workers:
            w.terminate()     

    return rsaResults, overlapResults


def collectRSASimpleResults(sitesInfo, resultsQueue, totalTasks, noSeqs, stream=None):
    stats = initSimpleRSAStats(sitesInfo)
    total = 0    
    while totalTasks:                
        result = resultsQueue.get()
        totalTasks -= 1                           
        if result is None:
            continue        
#         print(total, resultsQueue.qsize())      
        statsi = result        
        # update relevant statistics
        stats["seqsCutByAny"] += statsi["seqsCutByAny"]
        for site in sitesInfo.keys():
            stats["siteHitsCount"][site] += statsi["siteHitsCount"][site]
            stats["siteHitSeqsCount"][site] += statsi["siteHitSeqsCount"][site]
            stats["siteHitsSeqsIDs"][site] += statsi["siteHitsSeqsIDs"][site]
            # print(statsi["siteHitsSeqsIDs"][site])
        total += statsi["total"]
        if total % 50000 == 0:
            printto(stream, '\t%d/%d records have been collected ... ' % (total, noSeqs))
            sys.stdout.flush()
    printto(stream, '\t%d/%d sequences have been collected ... ' % (total, noSeqs))
    stats["total"] = noSeqs
    for site in sitesInfo.keys():
        stats["siteHitsSeqsIDs"][site] = set(stats["siteHitsSeqsIDs"][site])
    return stats
# for id in iter:
#         record = records[id]
#         try:              
#             seq = str(record.seq)
#             seqRC = str(Seq(seq).reverse_complement())
#             cut = False
#             for site in stats["siteHitsCount"].keys():
#                 hits = findHits(seq, sitesInfo[site])                
#                 if len(hits) == 0:
#                     hits = findHits(seqRC, sitesInfo[site])                   
#                 if len(hits) > 0:
#                     stats["siteHitsCount"][site] += len(hits) 
#                     stats["siteHitSeqsCount"][site] += 1                     
#                     stats["siteHitsSeqsIDs"][site].add(record.id)   
#                     cut = True                 
#             if cut:
#                 stats["seqsCutByAny"] += 1
#             procSeqs += 1
#             if procSeqs % seqsPerFile == 0:
#                 print('%d/%d sequences have been searched ... ' % (procSeqs, stats["total"]))
#                 sys.stdout.flush()
# #                 break
#         except BaseException as e:                
#             print(e)
#             raise        
#     print('%d/%d sequences have been searched ... ' % (procSeqs, stats["total"]))
#     records.close()


def loadRestrictionSites(sitesFile, stream=None):
    """
    given a whitespace separated file containing 2 columns, return a dictionary of restriction enzyme names to
    a regex translated sequence. Ignores all lines that starts with "#"

    :param sitesFile:
    :param stream:
    :return:
    """
    with open(sitesFile) as fp:
        sites = {}
        for line in fp:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    enzyme, seq = line.split()
                    if enzyme in sites:
                        printto(stream, enzyme + " is duplicated, the older enzyme sequence {} ".format(sites[enzyme]) +
                                                 "will be overridden.", LEVEL.WARN)
                    sites[enzyme] = replaceIUPACLetters(str(seq).upper().strip())
                except Exception as e:
                    printto(stream, "Offending line: {}, {}".format(line, line.split()), LEVEL.EXCEPT)
                    raise e

    printto(stream, "Restricting sites have been loaded")
    return sites


def replaceIUPACLetters(iupacSeq):
    """
    translates IUPAC letters to regex ACGT letters
    :param iupacSeq: string of IUPAC sequence
    :return: equivalent IUPAC sequence in a ACGT regex string
    """
    iupac = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': '[AG]',
        'Y': '[CT]',
        'S': '[GC]',
        'W': '[AT]',
        'K': '[GT]',
        'M': '[AC]',
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'V': '[ACG]',
        'N': '.'
    }
    tcgaSeq = ''
    for s in iupacSeq.upper():
        if s not in iupac:
            tcgaSeq += s
        else:
            tcgaSeq += iupac[s]
    return tcgaSeq


def findHitsRegion(qsRec, hitStarts):
    """
    return framework / cdr region where hitStarts is located at

    :param qsRec: dict
                row of cloneAnnot object (dict-like)

    :param hitStarts: list of ints
                indices of where the match starts

    :return: dict
                for each hit in a region, save the region as a key with value 1 (even if multiple hits are
                in the same region)
    """
    vhStart = qsRec['vqstart'] - qsRec['vstart']
    regions = {}
    for s in hitStarts:
        if ((qsRec['fr1.start'] - qsRec['vstart'] - vhStart) <= s <= (qsRec['fr1.end'] - vhStart)):
            regions['fr1'] = 1
        elif ((qsRec['cdr1.start'] - vhStart) <= s <= (qsRec['cdr1.end'] - vhStart)):
            regions['cdr1'] = 1
        elif ((qsRec['fr2.start'] - vhStart) <= s <= (qsRec['fr2.end'] - vhStart)):
            regions['fr2'] = 1
        elif ((qsRec['cdr2.start'] - vhStart) <= s <= (qsRec['cdr2.end'] - vhStart)):
            regions['cdr2'] = 1
        elif ((qsRec['fr3.start'] - vhStart) <= s <= (qsRec['fr3.end'] - vhStart)):
            regions['fr3'] = 1
        elif ((qsRec['cdr3.start'] - vhStart) <= s <= (qsRec['cdr3.end'] - vhStart)):
            regions['cdr3'] = 1
        elif ((not isnan(qsRec['fr4.end'])) and ((qsRec['fr4.start'] - vhStart) <= s <= (qsRec['fr4.end'] - vhStart))):
            regions['fr4'] = 1    
        else:
            raise Exception("Expected index {} to be within one of the FR/CDR regions. Record = {} vhStart = {}"
                            .format(s, qsRec, vhStart))
    return regions


def findHits(seq, site):
    seq = seq.upper()
    site = site.replace('/', '')  
    return [match.start() for match in re.finditer('(?=(%s))' % site, seq)]
#     return len(re.findall(site, seq))
