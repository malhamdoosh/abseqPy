'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

from numpy import isnan, nan
import re
from Bio import SeqIO
from multiprocessing import Queue, Manager
import sys
from math import ceil
from config import MEM_GB
from RestrictionSitesScanner import RestrictionSitesScanner
from productivityAuxiliary import ProcCounter
from pandas.core.frame import DataFrame
from IgRepAuxiliary.SeqUtils import readSeqFileIntoDict


def initSimpleRSAStats(sites):
    stats = {} # replace all with stats dict to be passed to a reporting function
    stats["siteHitsCount"] = {}
    stats["siteHitSeqsCount"] = {}
    stats["siteHitsSeqsIDs"] = {}
    stats["seqsCutByAny"] = 0
    stats["total"] = 0
    for site in sites.keys():
        stats["siteHitsCount"][site] = 0
        stats["siteHitSeqsCount"][site] = 0
        stats["siteHitsSeqsIDs"][site] = []
    return stats
        

def postProcessRSA(stats, sitesInfo):    
    rsaResults = []
    sites = sorted(stats["siteHitSeqsCount"], key=stats["siteHitSeqsCount"].get)    
    for site in sites:
        rsaResults.append([site, sitesInfo[site],
              stats["siteHitsCount"][site],
            stats["siteHitsCount"][site] * 100.0 / sum(stats["siteHitsCount"].values()),
            stats["siteHitSeqsCount"][site],
            stats["siteHitSeqsCount"][site] * 100.0 / stats["total"]])
    rsaResults.append(["Cut by any", nan, nan, nan,
                       stats["seqsCutByAny"], 
                       stats["seqsCutByAny"] * 100.0 / stats["total"]])
    rsaResults.append(["Total", nan, nan, nan,
                       stats["total"], 
                       100])
    rsaResults = DataFrame(rsaResults, 
               columns=["Enzyme", "Restriction Site" , "No.Hits" , "Percentage of Hits (%)" ,
                        "No.Molecules", "Percentage of Molecules (%)"])
    overlapResults = {}
    overlapResults["order1"] = stats["siteHitsSeqsIDs"]
    if len(stats["siteHitsSeqsIDs"]) >= 3:
        overlapResults["order2"] = calcRSAOverlapOrder2(stats["siteHitsSeqsIDs"], sites)
    return rsaResults, overlapResults

def calcRSAOverlapOrder2(order1, sites):
    print("The 2nd order overlapping matrix is being calculated using Jaccard Index ... ")
    overlap = []
    for site1 in sites:
        overlap.append([])
        for site2 in sites:
            inter = len(order1[site1].intersection(order1[site2])) * 1.0
            uni = len(order1[site1].union(order1[site2]))
            if (uni != 0):
                overlap[-1].append(inter / uni)
            else:
                overlap[-1].append(1)
    overlap = DataFrame(overlap, columns=sites, index=sites)    
    #overlap = linkage(overlap)  
    return overlap



def scanRestrictionSitesSimple(name, readFile, format, 
       cloneAnnot, sitesFile, threads):
    sitesInfo = loadRestrictionSites(sitesFile)
    seqsPerWorker = len(sitesInfo)
    workers = []   
    try:
        m = Manager()        
        records = m.dict()
        readSeqFileIntoDict(readFile, format, records)                
        queryIds = cloneAnnot.index
        noSeqs = len(queryIds)
        print("{0:,} restriction sites are being scanned for {1:,} sequences ...".format(
                                                         len(sitesInfo), noSeqs))
        ### Parallel implementation of the refinement        
        totalTasks = int(ceil(noSeqs * 1.0 / seqsPerWorker)) 
        tasks = Queue()      
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs, desc="sequences") 
        if (threads > totalTasks):
            threads = totalTasks     
        if MEM_GB < 16:
            threads = 2  
        # Initialize workers
        for i in range(threads):
            w = RestrictionSitesScanner(records, cloneAnnot, procCounter, 
                                        sitesInfo.copy(), simpleScan = True)
            w.tasksQueue = tasks
            w.exitQueue = exitQueue  
            w.resultsQueue = resultsQueue    
            workers.append(w)
            w.start()
        if (totalTasks > 1): 
            for i in range(totalTasks):
                #t = time.time()
                ids = queryIds[i * seqsPerWorker : (i+1) * seqsPerWorker] 
                tasks.put(ids)      
                #print("Preparing a job took: {0:f}".format(time.time() - t))          
        else:                
            recs = map(lambda x: records[x], queryIds)            
            tasks.put(recs)
        # Add a poison pill for each worker
        for i in range(threads + 10):
            tasks.put(None)       
        # Wait all process workers to terminate                
        i = 0 
        while i < threads:    
            m = exitQueue.get()
            if m == "exit":
                i += 1
        print("All workers have completed their tasks successfully.")
        # Collect results
        print("Results are being collated from all workers ...")  
        # invoking the result collection method   
        stats = collectRSASimpleResults(sitesInfo, resultsQueue, totalTasks, 
                   noSeqs)
        (rsaResults, overlapResults) = postProcessRSA(stats, sitesInfo)
        ### End of parallel implementation                  
        sys.stdout.flush()          

        print("Results were collated successfully.")
    except Exception as e:
        print("Something went wrong during the RSA scanning process!")
        raise e
    finally:        
        for w in workers:
            w.terminate()     
        #records.close()
    return (rsaResults, overlapResults)

def collectRSASimpleResults(sitesInfo, resultsQueue, totalTasks, 
                   noSeqs):
    stats = initSimpleRSAStats(sitesInfo)
    total = 0    
    while totalTasks:                
        result = resultsQueue.get()
        totalTasks -= 1                           
        if (result is None):
            continue        
#         print(total, resultsQueue.qsize())      
        statsi = result        
        # update relevant statistics
        stats["seqsCutByAny"] += statsi["seqsCutByAny"]
        for site in sitesInfo.keys():
            stats["siteHitsCount"][site] += statsi["siteHitsCount"][site]
            stats["siteHitSeqsCount"][site] += statsi["siteHitSeqsCount"][site]
            stats["siteHitsSeqsIDs"][site] +=  statsi["siteHitsSeqsIDs"][site]
            #print(statsi["siteHitsSeqsIDs"][site])
        total += statsi["total"]
        if (total % 50000 == 0):
            print('\t%d/%d records have been collected ... ' % (total, noSeqs))
            sys.stdout.flush()
    print('\t%d/%d sequences have been collected ... ' % (total, noSeqs))
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
    
def loadRestrictionSites(sitesFile):
    f = open(sitesFile)
    lines = f.readlines()
    sites = {}
    for line in lines:
        line = line.strip()
        if line and not line.startswith("#"):
            try:
                fields = line.split()
                if (sites.get(fields[0], None) is not None):
                    print(fields[0] + " is duplicated.")
                site = fields[1].upper()
                site = replaceIUPACLetters(site)
                site = site.replace('N', '.').replace('(', '[').replace(')', ']')                
                sites[fields[0]] = site
            except:
                print line, fields
                raise
    print("Restricting sites have been loaded")
    return sites


iupac = {
        'A':'A',
        'C': 'C',
        'G':'G',
        'T':'T',
        'R': '(AG)',
        'Y': '(CT)',
        'S': '(GC)',
        'W': '(AT)',
        'K': '(GT)',
        'M': '(AC)',
        'B': '(CGT)',
        'D': '(AGT)',
        'H': '(ACT)',
        'V': '(ACG)',
        'N':'N'         
        }
def replaceIUPACLetters(iupacSeq):
    tcgaSeq = ''
    iupacLetters = ''.join(iupac.keys())
    for s in iupacSeq.upper():
        if s not in iupacLetters:
            tcgaSeq += s
        else:
            tcgaSeq += iupac[s]
    return tcgaSeq
    
'''
    Used for restriction sites search 
'''
def findHitsRegion(cdrRec, hitStarts):
    vhStart = cdrRec['vqstart'] - cdrRec['vstart']
    regions = {}
    for s in hitStarts:
        if (s >= cdrRec['fr1.start']-cdrRec['vstart']  - vhStart and s <= cdrRec['fr1.end'] - vhStart):
            regions['fr1'] = 1
        elif (s >= cdrRec['cdr1.start'] - vhStart and s <= cdrRec['cdr1.end'] - vhStart):
            regions['cdr1'] = 1
        elif (s >= cdrRec['fr2.start'] - vhStart and s <= cdrRec['fr2.end'] - vhStart):
            regions['fr2'] = 1
        elif (s >= cdrRec['cdr2.start'] - vhStart and s <= cdrRec['cdr2.end'] - vhStart):
            regions['cdr2'] = 1
        elif (s >= cdrRec['fr3.start'] - vhStart and s <= cdrRec['fr3.end'] - vhStart):
            regions['fr3'] = 1
        elif (s >= cdrRec['cdr3.start'] - vhStart and s <= cdrRec['cdr3.end'] - vhStart):
            regions['cdr3'] = 1
        elif (not isnan(cdrRec['fr4.end']) and s >= cdrRec['fr4.start'] - vhStart and s <= cdrRec['fr4.end'] - vhStart):
            regions['fr4'] = 1    
        else:
            print(hitStarts, vhStart, cdrRec)
            raise
    return regions
                
def findHits(seq, site):
    seq = seq.upper()
    site = site.replace('/', '')  
    return [match.start() for match in re.finditer('(?=(%s))' %(site), seq)]
#     return len(re.findall(site, seq))




  