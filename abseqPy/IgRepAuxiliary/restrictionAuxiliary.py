'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import division
import re

from collections import defaultdict, Counter

from numpy import isnan, nan
from multiprocessing import Queue, Manager
from math import ceil
from pandas.core.frame import DataFrame

from abseqPy.IgRepAuxiliary.seqUtils import readSeqFileIntoDict
from abseqPy.IgRepAuxiliary.RestrictionSitesScanner import RestrictionSitesScanner
from abseqPy.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseqPy.logger import printto, LEVEL


def initRSAStats(simple):
    stats = {
        "siteHitsCount": defaultdict(int),
        "siteHitSeqsCount": defaultdict(int),
        "siteHitsSeqsIDs": defaultdict(set),
        "seqsCutByAny": 0,
        "total": 0,
    }

    # additional fields for detailed RSA
    if not simple:
        stats["siteHitsSeqsIGV"] = defaultdict(set)
        stats["hitRegion"] = defaultdict(Counter)
        stats["siteHitSeqsGermline"] = defaultdict(list)

    return stats


def calcRSAOverlapOrder2(order1, sites, stream=None):
    """
    returns a n by n matrix where n is len(sites) of jaccard index

    :param order1: dictionary of sets of ids
    :param sites: collection of enzymes
    :param stream: logging stream
    :return: n by n dataframe that has the form of a named matrix:

         enz1 enz2 enz3
    enz1    1  0.3  0.4
    enz2  0.3    1  0.5
    enz3  0.5  0.5    1
    """
    printto(stream, "The 2nd order overlapping matrix is being calculated using Jaccard Index ... ")
    overlap = []
    for site1 in sites:
        overlap.append([])
        for site2 in sites:
            inter = len(order1[site1].intersection(order1[site2]))
            uni = len(order1[site1].union(order1[site2]))
            if uni != 0:
                overlap[-1].append(inter / uni)
            else:
                overlap[-1].append(1)
    overlap = DataFrame(overlap, columns=sites, index=sites)
    # overlap = linkage(overlap)
    return overlap


def scanRestrictionSites(name, readFile, cloneAnnot, sitesFile, threads, simple=True, stream=None):
    """
    :param name: string
            analysis name

    :param readFile: string
            raw FASTQ/FASTA file

    :param cloneAnnot: dataframe
            IgRepertoire.cloneAnnot dataframe, depending on what the argument to 'simple' is, will require at least
            ['vqstart', 'vstart', 'fr4.end'] and 'vgene' if simple=False columns defined in the dataframe

    :param sitesFile: string
            path to restriction sites enzyme whitespace separated file. Example:
            enzyme1\tACGYTARRB
            enzyme2\tYTBABBAATG
            ...

    :param threads: int

    :param simple: bool
            simple or detailed analysis

    :param stream: logging stream

    :return: 2-tuple:
    (
        dataframe with columns: "Enzyme", "Restriction Site", "No.Hits", "Percentage of Hits (%)",
                                "No.Molecules", "Percentage of Molecules (%)"
                                where
                                    1. No.Hits are the total number of found hits for an enzyme (one molecule may have
                                        multiple enzyme hits)
                                    2. No.Molecules are the total number of molecules that the enzyme matched against.
                                       (if a molecule has multiple hotspots, only one is counted)
        dictionary with optional keys:
            {
                "order1" : {'enzyme1': {'seq_id1', 'seq_id2', 'seq_id3', ...}, 'enzyme2': {'seq_id5', ...} , ... },
                "order2" : Dataframe of n^2(all-vs-all) rows where each row is a jaccard index of the ids that each
                           pairwise comparison of the enzyme yields. This dataframe has an index column and header
                           that is identical (i.e. a "named matrix") - see calcRSAOverlapOrder2's return value
            }
            "order1" is always there, "order2" only appears if the number of enzymes is at least 3(len(sitesInfo)) >= 3)
    )
    """
    sitesInfo = loadRestrictionSites(sitesFile, stream=stream)
    seqsPerWorker = len(sitesInfo)
    workers = []
    try:
        m = Manager()
        records = m.dict()
        readSeqFileIntoDict(readFile, records, stream=stream)
        queryIds = cloneAnnot.index
        noSeqs = len(queryIds)
        printto(stream, "{:,} restriction sites are being scanned for {:,} sequences ..."
                .format(len(sitesInfo), noSeqs))
        totalTasks = int(ceil(noSeqs / seqsPerWorker))
        tasks = Queue()
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs, desc="sequences", stream=stream)
        threads = min(threads, totalTasks)

        # Initialize workers
        for _ in range(threads):
            # pass in a minified cloneAnnot, the bare minimum of what's required in the runSimple() / runDetailed()
            # method
            w = RestrictionSitesScanner(records, cloneAnnot, procCounter,
                                        sitesInfo.copy(), simpleScan=simple, stream=stream)
            w.tasksQueue = tasks
            w.exitQueue = exitQueue
            w.resultsQueue = resultsQueue
            workers.append(w)
            w.start()

        assert totalTasks > 0

        for i in range(totalTasks):
            tasks.put(queryIds[i * seqsPerWorker:(i + 1) * seqsPerWorker])

        # Add a poison pill for each worker
        for _ in range(threads + 10):
            tasks.put(None)

        # Wait all process workers to terminate
        i = 0
        while i < threads:
            i += (exitQueue.get() == "exit")

        # Collect results
        printto(stream, "All workers have completed their tasks successfully.")
        printto(stream, "Results are being collated from all workers ...")
        stats = collectRSAResults(sitesInfo, resultsQueue, totalTasks, noSeqs, simple=simple, stream=stream)
        (rsaResults, overlapResults) = postProcessRSA(stats, sitesInfo, simple=simple, stream=stream)
        printto(stream, "Results were collated successfully.")

    except Exception as e:
        printto(stream, "Something went wrong during the RSA scanning process, error: {}".format(str(e)), LEVEL.EXCEPT)
        raise e
    finally:
        for w in workers:
            w.terminate()

    return rsaResults, overlapResults


def collectRSAResults(sitesInfo, resultsQueue, totalTasks, noSeqs, simple=True, stream=None):
    stats = initRSAStats(simple=simple)
    total = 0
    while totalTasks:
        statsi = resultsQueue.get()
        if statsi is None:
            continue
        totalTasks -= 1

        # -------- update relevant statistics -------  #

        # 1. total number of sequences that are cut by any sites at all (i.e. number of sequences that are cut by
        #    *at least* one site)
        stats["seqsCutByAny"] += statsi["seqsCutByAny"]
        for site in sitesInfo.keys():
            # 2. total number of "possible hits" of this 'site' on all sequences (note, multi-hits are counted)
            stats["siteHitsCount"][site] += statsi["siteHitsCount"][site]
            # 3. total number of "hits" of this 'site' on all sequences (note, multi-hits on one sequence are still
            #    counted as one, not multi) - this is a "duplicate" field of siteHitsSeqsIDs, we could've taken the
            #    length of sitHitsSeqsIDs, it would be equal to this. This is left here from legacy code.
            stats["siteHitSeqsCount"][site] += statsi["siteHitSeqsCount"][site]
            # 4. the ids of which this site has at least one match, this length of this value should be equal to
            #    siteHitSeqsCount
            stats['siteHitsSeqsIDs'][site] = stats["siteHitsSeqsIDs"][site].union(statsi["siteHitsSeqsIDs"][site])

            if not simple:
                # these keys are only available for detailed RS analysis

                # 5. collect the total number of region where a match with this site has been registered
                # Counter object
                stats['hitRegion'][site] += statsi['hitRegion'][site]

                # 6. collect all the germline sequences that were recorded during a match with this site
                # list object
                stats['siteHitSeqsGermline'][site] += statsi['siteHitSeqsGermline'][site]

                # 7. collect all the IGV sequences that were recorded during a match with this site
                # set object
                stats['siteHitsSeqsIGV'][site] = stats['siteHitsSeqsIGV'][site].union(statsi['siteHitsSeqsIGV'][site])

        total += statsi["total"]
        if total % 50000 == 0:
            printto(stream, '\t%d/%d records have been collected ... ' % (total, noSeqs))

    printto(stream, '\t%d/%d sequences have been collected ... ' % (total, noSeqs))

    assert total == noSeqs
    stats["total"] = noSeqs
    return stats


def postProcessRSA(stats, sitesInfo, simple=True, stream=None):
    """
    returns a processed RSA result tuple. see return value for more information

    :param stats: dictionary of stats. see collectRSASimpleResults for the exact format
    :param sitesInfo: dictionary of enzymes mapped to their compiled regex
    :param simple: bool. simple or detailed RSA
    :param stream: logging stream
    :return: 2-tuple:
    (
        dataframe with columns: "Enzyme", "Restriction Site", "No.Hits", "Percentage of Hits (%)",
                                "No.Molecules", "Percentage of Molecules (%)",
                                <'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4', 'V Germlines'>
                                where
                                    1. No.Hits are the total number of found hits for an enzyme (one molecule may have
                                        multiple enzyme hits)
                                    2. No.Molecules are the total number of molecules that the enzyme matched against.
                                       (if a molecule has multiple hotspots, only one is counted)
                                    3. Columns in angle brackets <> are only present if simple = False (i.e. detailed
                                       RSA)
        dictionary with optional keys:
            {
                "order1" : {'enzyme1': {'seq_id1', 'seq_id2', 'seq_id3', ...}, 'enzyme2': {'seq_id5', ...} , ... },
                "order2" : Dataframe of n^2(all-vs-all) rows where each row is a jaccard index of the ids that each
                           pairwise comparison of the enzyme yields
            }
            "order1" is always there, "order2" only appears if the number of enzymes is at least 3(len(sitesInfo)) >= 3)
    )
    """
    rsaResults = []
    extraColumns = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4', 'V Germlines']
    # sort site names by their highest sequence hit (deduplicated) count (i.e. non-multi-hit count)
    sites = sorted(stats["siteHitSeqsCount"], key=stats["siteHitSeqsCount"].get)
    for site in sites:
        # site name, site regex pattern, total number of hits (incl multi-hit), normalized total multi-hit, \
        #       total number molecules that are hit by this site, normalized tot. no.mol that are hit by site (dedup)
        rowData = [site, sitesInfo[site].pattern, stats["siteHitsCount"][site],
                   (stats["siteHitsCount"][site] / sum(stats["siteHitsCount"].values())) * 100,
                   stats["siteHitSeqsCount"][site],
                   (stats["siteHitSeqsCount"][site] / stats["total"]) * 100]
        if not simple:
            # add FR1, CDR1, ... , FR4 region stats, add IGV stats
            rowData += [stats["hitRegion"][site][region] for region in extraColumns[:-1]] + \
                       ['|'.join(stats["siteHitsSeqsIGV"][site])]
            # # todo(ASK) 1: paired with ask 1. Find the other TODO to see what this is about
            # seqs = []
            # for strand, seq in stats['siteHitSeqsGermline'][site]:
            #     seqs.append(SeqRecord(Seq(seq), id=('seq' + str(len(seqs)) + strand)))
            # SeqIO.write(seqs, siteHitsFile.replace(".csv", "_germline" + site + ".fasta"), "fasta")

        rsaResults.append(rowData)

    # total number of sequences that are cut by at least one site (and normalized)
    rsaResults.append(
        ["Cut by any", nan, nan, nan, stats["seqsCutByAny"], stats["seqsCutByAny"] / stats["total"]] +
        ([] if simple else [nan] * len(extraColumns))  # we have extra columns in detailed RSA
    )

    # total number of sequences (and normalized, 100%, doh!)
    rsaResults.append(["Total", nan, nan, nan, stats["total"], 100] +
                      ([] if simple else [nan] * len(extraColumns))  # we have extra columns in detailed RSA
                      )

    rsaResults = DataFrame(rsaResults,
                           columns=(["Enzyme", "Restriction Site", "No.Hits",
                                    "Percentage of Hits (%)", "No.Molecules", "Percentage of Molecules (%)"] +
                                    ([] if simple else extraColumns))
                           )

    overlapResults = {"order1": stats["siteHitsSeqsIDs"]}

    # if there are at least 3 sites, calculate overlap order2
    if len(stats["siteHitsSeqsIDs"]) >= 3:
        overlapResults["order2"] = calcRSAOverlapOrder2(stats["siteHitsSeqsIDs"], sites, stream=stream)

    return rsaResults, overlapResults


def loadRestrictionSites(sitesFile, stream=None):
    """
    given a whitespace separated file containing 2 columns, return a dictionary of restriction enzyme names to
    a regex translated sequence. Ignores all lines that starts with "#"

    :param sitesFile: file with 2 cols, enzyme <ws> seq. Any line that *starts* with # will be ignored
    :param stream: logging stream
    :return: dictionary of enzyme to precompiled regex mapping, for example:
    {
        "ENZYME1": re.compile("AC[GT]..A")     # assuming ENZYME1's IUPAC sequence was "ACKNNA"
    }
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
                    sites[enzyme] = re.compile(replaceIUPACLetters(str(seq).upper().strip()))
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
    return framework / cdr region where hitStarts is located at. Raises exception of any one of hitStart's elements
    don't fall between fr1.start and fr4.end

    :param qsRec: dict
                row of cloneAnnot object (dict-like)

    :param hitStarts: list of ints
                indices of where the match starts across the sequence represented by qsRec (0-based index)

    :return: dict
                for each hit in a region, save the region as a key with value 1 (even if multiple hits are
                in the same region)
    """
    vhStart = max(0, qsRec['vqstart'] - qsRec['vstart'])
    regions = {}
    for s in hitStarts:
        # offset index by the number of nucleotides that were trimmed when matching for restriction sites enzyme
        # +1 since s is 0-based index
        s += (vhStart + 1)

        if qsRec['fr1.start'] <= s <= qsRec['fr1.end']:
            regions['fr1'] = 1
        elif qsRec['cdr1.start'] <= s <= qsRec['cdr1.end']:
            regions['cdr1'] = 1
        elif qsRec['fr2.start'] <= s <= qsRec['fr2.end']:
            regions['fr2'] = 1
        elif qsRec['cdr2.start'] <= s <= qsRec['cdr2.end']:
            regions['cdr2'] = 1
        elif qsRec['fr3.start'] <= s <= qsRec['fr3.end']:
            regions['fr3'] = 1
        elif qsRec['cdr3.start'] <= s <= qsRec['cdr3.end']:
            regions['cdr3'] = 1
        elif (not isnan(qsRec['fr4.end'])) and (qsRec['fr4.start'] <= s <= qsRec['fr4.end']):
            regions['fr4'] = 1
        else:
            raise Exception("Expected index {} to be within one of the FR/CDR regions. Record = {} vhStart = {}"
                            .format(s, qsRec, vhStart))
    return regions


def findHits(seq, site):
    """
    returns non overlapping matching indices of "site" in "seq"

    :param seq: nucleotide string
    :param site: compiled site, for example, re.compile("AC[GT].T")
    :return: a list of start indices, where "seq" matched "site"
    """
    seq = seq.upper()
    return [match.start() for match in site.finditer(seq)]
