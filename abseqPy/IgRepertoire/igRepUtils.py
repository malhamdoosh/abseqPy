'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import print_function

import gzip
import shutil
import glob
import re
import os
import sys

from os.path import exists
from Bio import SeqIO, AlignIO
from subprocess import CalledProcessError
from pandas.core.frame import DataFrame
from numpy import isnan, nan
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter, defaultdict
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat import MatrixInfo as matlist

from abseqPy.config import CLUSTALOMEGA, IGBLASTN, IGBLASTP, LEEHOM, PEAR, FLASH
from abseqPy.logger import printto, LEVEL
from abseqPy.utilities import hasLargeMem, ShortOpts, quote


def detectFileFormat(fname, noRaise=False):
    """
    detects if the filename ends with fastq or fasta extensions (it can be zipped)
    :param fname: filename for which the extension should be identified (fname can be zipped)
    :return: "fastq" or "fasta" or None depending on the extensions. This output is compatible with the expected
    input argument to the 'format' parameter in SeqIO.* function variants
    """

    class FileFormatNotSupported(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    if ".fastq" in fname or ".fq" in fname:
        return "fastq"
    if ".fasta" in fname or ".fa" in fname:
        return "fasta"
    if not noRaise:
        raise FileFormatNotSupported("Only FASTQ or FASTA (.fastq, .fq, .fasta, .fa) extensions are supported")
    return None


def inferSampleName(fname, merger, fastqc):
    """
    infers the sample name from a given file.
    EG: SRR1002_L001.fastq.gz   => SRR1002_L001
        Sample1_L001.fastq.gz   => Sample1_L001
        Sample1.fastq.gz        => Sample1

    :param fname: string
            filename to infer from (usually args.f1)

    :return: string
            inferred sample name from fname
    """
    f1name = os.path.basename(fname)
    # read is paired end ==> remove everything after _R1.fast*...
    if f1name.find("_R") != -1 and (merger or fastqc):
        ext = '_' + f1name.split("_")[-1]
    else:
        ext = f1name[f1name.find("."):]
    sampleName = f1name.replace(ext, "")
    # sampleName = f1name.split("_")[0] + "_"
    # sampleName += f1name.split("_")[-1].split(".")[0]
    return sampleName


# U flag = Universal ending flag (windows/dos/mac/linux  ... etc) (http://biopython.org/wiki/SeqIO)
def safeOpen(filename, mode="rU"):
    """
    given a filename, if it is ending with a gzipped extension, open it with gzip, else open normally
    :param filename: file to be opened
    :param mode: mode to open file in
    :return: file handle
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def gunzip(gzipFile):
    """
    Given a gzipped file, create a similar file that's uncompressed. If the file is not gzipped, do nothing.
    The naming scheme follows the original provided path to file, but with .gz suffix stripped
    The original gzipped file stays as a zipped file.

    NOTE: WILL NOT DECOMPRESS IF THE DECOMPRESSED FILE IS FOUND

    :param gzipFile: file(filename) to be unzipped
    :return: new filename of uncompressed file, or if file was originally not gzipped, return same name as argument
    """

    if not gzipFile.endswith(".gz"):
        return gzipFile

    newFileName = gzipFile.replace(".gz", "")
    # only decompress if the file is not found
    if not os.path.exists(newFileName):
        with gzip.open(gzipFile, 'rb') as f_in, open(newFileName, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return newFileName


def fastq2fasta(fastqFile, outputDir, stream=None):
    """
    Converts a fastq file into fasta file. Fastq can be compressed if it was provided as such
    :param fastqFile: (un)compressed fastq file. If compressed, will leave original compressed untouched
    :param outputDir: Where to produce the new fasta file
    :param stream: debugging stream
    :return: fasta filename
    """
    # FASTQ to FASTA
    # awk 'NR % 4 == 1 {print ">" $0 } NR % 4 == 2 {print $0}' my.fastq > my.fasta
    filename = os.path.basename(fastqFile)
    seqOut = os.path.join(outputDir, "seq")

    if not os.path.isdir(seqOut):
        os.makedirs(seqOut)

    # rename all fastq files to fasta, including gzipped files
    if filename.endswith(".gz"):
        filename = os.path.join(seqOut, filename.replace(filename.split('.')[-2] + ".gz", 'fasta'))
        fastqFile = gunzip(fastqFile)
    else:
        filename = os.path.join(seqOut, filename.replace(filename.split('.')[-1], 'fasta'))

    if exists(filename):
        printto(stream, "\tThe FASTA file was found!", LEVEL.WARN)
        return filename

    printto(stream, "\t" + os.path.basename(fastqFile) + " is being converted into FASTA ...")
    SeqIO.convert(fastqFile, 'fastq', filename, 'fasta')

    # not all systems have AWK by default (cough, windows)
    # command = ("awk 'NR % 4 == 1 {sub(\"@\", \"\", $0) ; print \">\" $0} NR % 4 == 2 "
    #            "{print $0}' " + fastqFile + " > " + filename
    #            )
    # os.system(command)
    return filename


def runIgblastn(blastInput, chain, threads=8,
                db='$IGBLASTDB', igdata="$IGDATA", domainSystem='imgt',
                outputDir="", species='human', stream=None):
    """
    Excerpt:
        The V domain can be delineated using either IMGT system (Lefranc et al 2003) or
        Kabat system (Kabat et al, 1991, Sequences of Proteins of Immunological Interest,
        National Institutes of Health Publication No. 91-3242, 5th ed., United States Department
        of Health and Human Services, Bethesda, MD).
        Domain annotation of the query sequence is based on pre-annotated domain information
        for the best matched germline hit.

        IMGT classification system is used to delineate the V domain

    this function simply runs the igblastn command and returns the output file name

    :param blastInput: path to input fasta file
    :param chain: chain type, one of hv, kv, lv, or klv
    :param threads: int
    :param db: path to the directory containing imgt_<species>_ig[hkl][vdj]. Will also read
    the environment variable $IGBLASTDB if not provided
    :param igdata: path to the directory containing optional and internal data used by IgBLAST
    This env variable is always defined otherwise igblast wouldn't work anyway
    :param domainSystem: string, one of imgt or kabat
    :param outputDir: output directory where the output file will be saved in
    :param species: human is the only species supported currently
    :param stream: logger stream object
    :return: output filename
    """
    if outputDir:
        head, tail = os.path.split(blastInput)
        blastOutput = os.path.join(outputDir, tail.replace('.' + tail.split('.')[-1], '.out'))
    else:
        blastOutput = blastInput.replace('.' + blastInput.split('.')[-1], '.out')

    if exists(blastOutput):
        printto(stream, "\tBlast results were found ... " + os.path.basename(blastOutput))
        return blastOutput

    printto(stream, '\tRunning igblast ... ' + os.path.basename(blastInput))

    command = buildIgBLASTCommand(igdata, db, chain, species, domainSystem,
                                  blastInput, blastOutput, threads, protein=False, stream=stream)
    try:
        command(stderr=sys.stderr, stdout=sys.stdout)
    except CalledProcessError as e:
        printto(stream, "Command {} failed with error code {}.\nDUMP: {}".format(command, e.returncode, e.output),
                LEVEL.CRIT)
        raise e
    return blastOutput


def runIgblastp(blastInput, chain, threads=8, db='$IGBLASTDB', igdata='$IGDATA', domainSystem='imgt',
                outputDir="", species='human', stream=None):
    """
    Excerpt:
        The V domain can be delineated using either IMGT system (Lefranc et al 2003) or
        Kabat system (Kabat et al, 1991, Sequences of Proteins of Immunological Interest,
        National Institutes of Health Publication No. 91-3242, 5th ed., United States Department
        of Health and Human Services, Bethesda, MD).
        Domain annotation of the query sequence is based on pre-annotated domain information
        for the best matched germline hit.

        IMGT classification system is used to delineate the V domain

    this function simply runs the igblastn command and returns the output file name

    :param blastInput: path to input fasta file
    :param chain: chain type, one of hv, kv, or lv
    :param threads: int
    :param db: path to the directory containing imgt_<species>_ig[hkl][vdj]. Will also read
    the environment variable $IGBLASTDB if not provided
    :param igdata: path to the directory containing optional and internal data used by IgBLAST
    This env variable is always defined otherwise igblast wouldn't work anyway
    :param domainSystem: string, one of imgt or kabat
    :param outputDir: output directory where the output file will be saved in
    :param species: human is the only species supported currently
    :param stream: logger stream object
    :return: output filename
    """
    blastOutput = os.path.join(outputDir, blastInput.replace('.' + blastInput.split('.')[-1], '.out'))

    if exists(blastOutput):
        printto(stream, "\tBlast results were found ... " + os.path.basename(blastOutput))
        return blastOutput

    printto(stream, '\tRunning igblast ... ' + os.path.basename(blastInput))

    command = buildIgBLASTCommand(igdata, db, chain, species, domainSystem, blastInput, blastOutput,
                                  threads, protein=True, vOnly=True, stream=stream)
    try:
        command(stderr=sys.stderr, stdout=sys.stdout)
    except CalledProcessError as e:
        printto(stream, "Command {} failed with error code {}.\nDUMP: {}".format(command, e.returncode, e.output),
                LEVEL.CRIT)
        raise e
    return blastOutput


def writeClonoTypesToFile(clonoTypes, filename, top=100, overRepresented=True, stream=None):
    if exists(filename):
        printto(stream, "\tThe clonotype file " + os.path.basename(filename) + " was found!", LEVEL.WARN)
        return

    total = sum(clonoTypes.values()) * 1.0
    dic = defaultdict(list)
    t = 0
    for k in sorted(clonoTypes, key=clonoTypes.get, reverse=overRepresented):
        dic['Clonotype'].append(str(k))
        dic['Count'].append(clonoTypes[k])
        dic['Percentage (%)'].append(clonoTypes[k] / total * 100)
        t += 1
        if t >= top:
            break

    df = DataFrame(dic)
    # fixed format (fast read/write) sacrificing search
    # (should change to table format(t) if search is needed for clonotype clustering/comparison)
    # df.to_hdf(filename, "clonotype", mode="w", format="f")
    df.to_csv(filename + ".gz", mode="w", compression="gzip")
    printto(stream, "\tA clonotype file has been written to " + os.path.basename(filename))


def writeCountsToFile(dist, filename):
    # This function prints the distribution counts into a text file
    with open(filename, 'w') as out:
        out.write('Germline group,Count,Percentage (%)\n')
        total = sum(dist.values()) * 1.0
        for k in sorted(dist, key=dist.get, reverse=True):
            out.write(str(k) + ',' + str(dist[k]) + ',' + ("%.2f" % (dist[k] / total * 100)) + '\n')
        out.write('TOTAL, ' + str(total) + ', 100 ')
    print("A text file has been created ... " + os.path.basename(filename))


def findBestAlignment(seq, query, dna=False, offset=0, show=False):
    if not dna:
        alignments = align.localds(seq.replace('*', 'X'), query, matlist.blosum62, -100, -100)
    else:
        alignments = align.localms(seq, query, 1, -2, -2, -2)

    #     print(seq, query, alignments)
    scores = [a[2] for a in alignments]
    if len(scores) == 0:
        return -1, -1, True
    best = scores.index(max(scores))
    if show:
        print(format_alignment(*alignments[best]))
        print(alignments[best])

    # FR4 start is where both sequence start to align with each other
    # including leading mismatches (these mismatches maybe due to mutations)
    #     0123456
    # eg: GGGGACGTACGTACGT
    #           ||||||||||
    #     ----CAGTACGTACGT
    # although alignment starts at pos 6, we still consider FR4 to start at pos 4
    start = extend5align(alignments[best]) + offset + 1     # 1-based start

    end = int(offset + alignments[best][-1])                # 1-based end

    gapped = False

    # subtract away non-existing '-'s from the seq because seq itself doesn't have these '-'s
    # eg: -GGGACGTACGTACGT
    #      |||||||||||||||
    #     GGGACAGTACGTACGT
    # should start at 1, not 2. because the leading '-' doesn't exist in the actual sequence!
    if '-' in alignments[best][0]:
        start -= alignments[best][0][:(alignments[best][-2] + 1)].count('-')
        end -= alignments[best][0][:(alignments[best][-1] + 1)].count('-')
        gapped = True
    return start, end, gapped  # 1-based


def extractProteinFrag(protein, start, end, offset=0, trimAtStop=False, stream=None):
    """
    Extract a protein fragment from a protein sequence based on DNA positions
    start and end are 1-based

    :param protein:
    :param start:
    :param end:
    :param offset:
    :param trimAtStop:
    :param stream:
    :return:
    """
    if isnan(start) or isnan(end):
        return ''
    if start != -1 and end != -1 and end - start < 1:
        return ''
    # start and end are 1-based positions
    start = (start - offset) if start != -1 else start
    end = (end - offset) if end != -1 else end
    try:
        if start != -1:
            # s = int(round((start  - 1.0 ) / 3))# 0-based
            s = int(((start - 1) / 3))  # 0-based
        else:
            s = 0
        if end != -1:
            # e = int(round( (end*1.0)  / 3)) # 1-based
            e = int(((end) / 3))  # 1-based
        else:
            e = len(protein)
        if (s + 1) < e:
            frag = protein[s:e]
        elif (s + 1) == e:
            frag = protein[s]
        else:
            return ''
        if trimAtStop and ('*' in frag):
            frag = frag[:frag.index('*')]
        return frag
    except:
        printto(stream, "ERROR at Extract Protein Fragment {} {} {}".format(protein, start, end), LEVEL.ERR)
        return None


def extractCDRsandFRsProtein(protein, qsRec, offset, stream=None):
    try:
        seqs = []
        newProtein = ""
        # Extract protein sequence of FR1             
        seqs.append(extractProteinFrag(protein, qsRec['fr1.start'], qsRec['fr1.end'], offset, stream=stream))
        # Extract protein sequence of CDR1
        seqs.append(extractProteinFrag(protein, qsRec['cdr1.start'], qsRec['cdr1.end'], offset, stream=stream))
        # Extract protein sequence of FR2
        seqs.append(extractProteinFrag(protein, qsRec['fr2.start'], qsRec['fr2.end'], offset, stream=stream))
        # Extract protein sequence of CDR2
        seqs.append(extractProteinFrag(protein, qsRec['cdr2.start'], qsRec['cdr2.end'], offset, stream=stream))
        # Extract protein sequence of FR3
        seqs.append(extractProteinFrag(protein, qsRec['fr3.start'], qsRec['fr3.end'], offset, stream=stream))
        # Extract protein sequence of CDR3 and FR4
        seqs.append(extractProteinFrag(protein, qsRec['cdr3.start'], qsRec['cdr3.end'], offset, stream=stream))
        seqs.append(extractProteinFrag(protein, qsRec['fr4.start'], qsRec['fr4.end'], offset, stream=stream))
        # check whether FR and CDR sequences were extracted correctly
        newProtein = ''.join(seqs)
        assert newProtein in protein
    except Exception as e:
        #         print("ERROR at partitioning the protein sequence: ")
        #         print (protein, newProtein, seqs, offset)
        #         raise e
        raise Exception("ERROR at partitioning the protein sequence: ")
    return newProtein, seqs


def extractCDRsandFRsDNA(dna, qsRec):
    try:
        seqs = []
        newDna = ""
        seqs.append(dna[int(qsRec['fr1.start'] - 1):int(qsRec['fr1.end'])])
        seqs.append(dna[int(qsRec['cdr1.start'] - 1):int(qsRec['cdr1.end'])])
        seqs.append(dna[int(qsRec['fr2.start'] - 1):int(qsRec['fr2.end'])])
        seqs.append(dna[int(qsRec['cdr2.start'] - 1):int(qsRec['cdr2.end'])])
        seqs.append(dna[int(qsRec['fr3.start'] - 1):int(qsRec['fr3.end'])])
        if (isnan(qsRec['cdr3.start']) or isnan(qsRec['cdr3.end']) or
                qsRec['cdr3.end'] - qsRec['cdr3.start'] < 0):
            seqs.append('')
        else:
            seqs.append(dna[int(qsRec['cdr3.start'] - 1):int(qsRec['cdr3.end'])])
        if isnan(qsRec['fr4.start']) or isnan(qsRec['fr4.end']):
            seqs.append('')
        else:
            seqs.append(dna[int(qsRec['fr4.start'] - 1):int(qsRec['fr4.end'])])
        newDna = ''.join(seqs)
        assert newDna in dna
    except Exception as e:
        #         print("ERROR at partitioning the nucleotide sequence: ")
        #         print (dna, newDna, seqs)
        #         raise e
        raise Exception("ERROR at partitioning the nucleotide sequence: ")
    return seqs


def mergeReads(readFile1, readFile2, threads=3, merger='leehom', outDir="./", stream=None):
    seqOut = os.path.join(outDir, "seq")

    if not os.path.isdir(seqOut):
        os.makedirs(seqOut)

    readFile = os.path.basename(readFile1)
    outputPrefix = os.path.join(seqOut, readFile.replace("_" + readFile.split('_')[-1], ''))

    if merger == 'pear':
        mergedFastq = outputPrefix + '.assembled.fastq'
        if not exists(mergedFastq):
            printto(stream, "{} and {} are being merged ...".format(os.path.basename(readFile1)
                                                                    , os.path.basename(readFile2)))
            pear = ShortOpts(exe=PEAR, f=quote(readFile1), r=quote(readFile2), o=quote(outputPrefix),
                             j=threads, v=15, n=350)
            # printto(stream, "Executing: " + str(pear))
            pear()
        else:
            printto(stream, "\tMerged reads file " + os.path.basename(mergedFastq) + ' was found!', LEVEL.WARN)
    elif merger == 'leehom':
        mergedFastq = outputPrefix + '.fq'
        if not exists(mergedFastq):
            printto(stream, "{} and {} are being merged ...".format(os.path.basename(readFile1)
                                                                    , os.path.basename(readFile2)))
            leehom = ShortOpts(LEEHOM, "--ancientdna", fq1=quote(readFile1), fq2=quote(readFile2),
                               fqo=quote(outputPrefix), t=threads)
            # printto(stream, "Executing: " + str(leehom))
            leehom()
            gunzip(mergedFastq + '.gz')
        else:
            printto(stream, "\tMerged reads file " + os.path.basename(mergedFastq) + ' was found!', LEVEL.WARN)
    elif merger == 'flash':
        mergedFastq = outputPrefix + '.extendedFrags.fastq'
        outputPrefix = os.path.basename(outputPrefix)
        if not exists(mergedFastq):
            printto(stream, "{} and {} are being merged ...".format(os.path.basename(readFile1)
                                                                    , os.path.basename(readFile2)))
            flash = ShortOpts(exe=FLASH, t=threads, o=quote(outputPrefix), r=300, f=450, s=50)\
                .append("{} {}".format(quote(readFile1), quote(readFile2)))
            # printto(stream, "Executing: " + str(flash))
            flash()
            for f in glob.glob("{}.*".format(outputPrefix)):
                shutil.move(f, seqOut)
        else:
            printto(stream, "\tMerged reads file " + os.path.basename(mergedFastq) + ' was found!', LEVEL.WARN)
    #     elif (merger == 'seqprep'):
    #         ### MERGE using SeqPrep
    #         mergedFastq = readFile1.replace(readFile1.split('_')[-1], 'merged.fastq.gz')
    #         unmerged1 = readFile1.replace('.fastq', '_unmerged.fastq.gz')
    #         unmerged2 = readFile2.replace('.fastq', '_unmerged.fastq.gz')
    #         aligns = readFile1.replace(readFile1.split('_')[-1], 'aligns.txt.gz')
    #         command = "SeqPrep -f %s -r %s -s %s -1 %s -2 %s -E %s"
    #         os.system(command % (readFile1, readFile2, mergedFastq,
    #                              unmerged1, unmerged2, aligns))
    #         os.system("gunzip " + mergedFastq)
    #         mergedFastq = mergedFastq.replace('.gz', '')
    #         ### END MERGE using SeqPrep
    else:
        raise Exception("Unknown short reads merger is selected")

    return os.path.abspath(mergedFastq)


def writeTableIntoFile(table, filename, stream=None):
    df = DataFrame(table)
    df.to_csv(filename, sep='\t', header=True, index=True)
    printto(stream, "Text file has been written to " + filename)


def writeListToFile(items, filename):
    out = open(filename, 'w')
    out.write("\n".join(items))
    out.close()


def compressSeqGeneLevel(seqDict):
    geneLevel = {}
    for ighv in seqDict.keys():
        gene = ighv.split('*')[0]
        if geneLevel.get(gene, None) is None:
            geneLevel[gene] = []
        geneLevel[gene] += seqDict[ighv]
    return geneLevel


def compressSeqFamilyLevel(seqDict):
    familyLevel = {}
    for ighv in seqDict.keys():
        fam = ighv.split('-')[0].split('/')[0]
        if familyLevel.get(fam, None) is None:
            familyLevel[fam] = []
        familyLevel[fam] += seqDict[ighv]
    return familyLevel


def compressCountsGeneLevel(countsDict):
    geneLevel = Counter()
    for k in countsDict.keys():
        ksub = k.split('*')[0]
        geneLevel[ksub] = geneLevel.get(ksub, 0) + countsDict[k]
    return geneLevel


def compressCountsFamilyLevel(countsDict):
    familyLevel = Counter()
    for k in countsDict.keys():
        ksub = k.split('-')[0].split('/')[0].rstrip('D')
        familyLevel[ksub] = familyLevel.get(ksub, 0) + countsDict[k]
    return familyLevel


def alignListOfSeqs(signals, outDir, threads, name, stream=None):
    """
    perform multiple sequence alignment using CLUSTAL

    :param signals:
    :param outDir:
    :param threads:
    :param name:
    :param stream:
    :return:
    """
    L = map(len, signals)

    printto(stream,
            "\t\t{} sequences are being aligned using CLUSTAL-OMEGA (L in [{}, {}])... ".format(len(L), min(L), max(L)))
    tempSeq = os.path.join(outDir, "csl_temp_seq_" + name + ".fasta")
    tempAlign = tempSeq.replace('.fasta', '.aln')

    seqs = []

    for i in range(len(signals)):
        seqs.append(SeqRecord(Seq(signals[i]), id='seq' + str(i)))
    SeqIO.write(seqs, tempSeq, 'fasta')

    clustal = ShortOpts(CLUSTALOMEGA, i=quote(tempSeq), o=quote(tempAlign))\
        .append("--threads={} --outfmt=clustal".format(threads))

    # printto(stream, "Executing: " + str(clustal))
    # throw away stderr and stdout
    clustal(stdout=None, stderr=None)

    alignment = AlignIO.read(tempAlign, 'clustal')
    alignedSeq = []
    for rec in alignment:
        alignedSeq.append(str(rec.seq))
    os.remove(tempSeq)
    os.remove(tempAlign)
    return alignedSeq


# source ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
matStr1 = ("   A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N,"
           "A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2,"
           "T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2,"
           "G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2,"
           "C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2,"
           "S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1,"
           "W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1,"
           "R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1,"
           "Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1,"
           "K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1,"
           "M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1,"
           "B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1,"
           "V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1,"
           "H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1,"
           "D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1,"
           "N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1"
           )
# custom: ambitious code matching == max score
matStr2 = ("   A   T   G   C,"
           "A   5  -4  -4  -4,"
           "T  -4   5  -4  -4,"
           "G  -4  -4   5  -4,"
           "C  -4  -4  -4   5,"
           "S  -4  -4   5   5,"
           "W   5   5  -4  -4,"
           "R   5  -4   5  -4,"
           "Y  -4   5  -4   5,"
           "K  -4   5   5  -4,"
           "M   5  -4  -4   5,"
           "B  -4  5  5  5,"
           "V  5  -4  5  5,"
           "H  5  5  -4  5,"
           "D  5  5  5  -4,"
           "N  5  5  5  5"
           )


def getIUPACSubMatrix():
    lines = matStr2.split(',')
    colHeads = lines[0].split()
    iupacSubMat = {}
    for line in lines[1:]:
        row = line.split()
        for i in range(len(row) - 1):
            iupacSubMat[(colHeads[i], row[0])] = float(row[i + 1])
    return iupacSubMat


subMatIUPAC = getIUPACSubMatrix()


def calMaxIUPACAlignScores(seqs):
    scores = []
    for seq in seqs:
        scores.append(0)
        for s in seq:
            rowscores = []
            for (r, c) in subMatIUPAC.keys():
                if (c == s):
                    rowscores.append(subMatIUPAC[(r, c)])
            scores[-1] += max(rowscores)
    return scores


def findBestMatchedPattern(seq, patterns, extend5end=False):
    """
    find the best matched pattern in a list of patterns
    and classify the type of the alignment (intact, indelled, mismatched, unknown)
    :param seq: nucleotide sequence
    :param patterns: zip iterator (or list) of (pattern_id, pattern_seq, pattern_max_IUPAC_score)
    :param extend5end: since this function uses Local edit distance, it will not favor mismatches and gaps earlier
    than the alignment. Use this flag to get the 'absolute beginning' of match
    :return: tuple of (pattern_id, mismatch_position, indel_position, start_pos (inclusive), end_pos (exclusive)).
    for example: (Oligo1H, 0, 0, 0, 15) means pattern id Oligo1H has the best match with 0 indel/mismatches and
    alignment starts from index 0 until 15: primer_seq[0:15]. If no alignment is ideal, returns (str(nan), 0, 0, -1, -1)

    Note: 0) mismatch_position and indel_position are 1-based index (i.e. starts from 1, not 0) - 0 means no indel/mis
          1) primer_id = 'nan'        => there was no suitable hit - mismatches and indel_pos will be left 0, but you
                                         should (obviously) not interpret that as mismatch at pos 0 or indel at pos 0
          2) mismatch_position = 0    => no mismatches
          3) indel_position = 0       => no indel_position
    """
    NO_MATCH = (str(nan), 0, 0, -1, -1)
    scores = []
    # align the sequence against all possible patterns
    for (id, pattern, maxScore) in patterns:
        alignments = align.localds(seq.upper(), pattern, subMatIUPAC, -5, -5)
        if len(alignments) > 1:
            localScores = [a[2] for a in alignments]
            alignment = alignments[localScores.index(max(localScores))]
        elif len(alignments) > 0:
            alignment = alignments[0]
        else:
            return NO_MATCH
        if alignment:
            alignLen = alignment[-1] - alignment[-2]
            scores.append((id, alignment))
            # if the sequence exactly matches one of the patterns (i.e. got the
            # max possible score from the matrix) ==> intact, return immediately
            if (alignment[2] == maxScore and
                    alignLen == len(pattern) and
                    '-' not in alignment[0] and
                    '-' not in alignment[1]):
                return scores[-1][0], 0, 0, alignment[-2], alignment[-1]
        else:
            scores.append((id, ('', '', 0)))

    # if no exact matching ==> find the best alignment (pattern)
    if len(scores) > 1:
        tmp = map(lambda x: x[1][2], scores)
        bestInd = tmp.index(max(tmp))
    elif len(scores) == 1:
        bestInd = 0
    else:
        return NO_MATCH

    best = list(scores[bestInd])
    best[1] = list(best[1])

    # best = [id, [seq, pattern, score, matchstart, matchend]]
    ID, ALIGNMENT = range(2)
    SEQ, PTN, SCORE, MSTART, MEND = range(5)

    if best[ALIGNMENT][SCORE] == 0:
        return NO_MATCH

    # classify the alignment type ==> insertion, deletion, mismatches

    # Find the position of Indel/Mismatch
    # remove starting indels
    if best[ALIGNMENT][PTN].startswith('-'):
        i = 0
        while best[ALIGNMENT][PTN][i] == '-':
            i += 1
        best[ALIGNMENT][SEQ] = best[ALIGNMENT][SEQ][i:]
        best[ALIGNMENT][PTN] = best[ALIGNMENT][PTN][i:]

    # TODO: revise algorithm
    # find the location of insertion or deletion 
    delPos = -1
    if '-' in best[ALIGNMENT][SEQ]:
        delPos = best[ALIGNMENT][SEQ].index('-')
    # if there is a gap at the beginning ==> happened because of insertion/deletion in the middle
    if '-' in best[ALIGNMENT][PTN] and best[ALIGNMENT][PTN].index('-') > delPos \
            and best[ALIGNMENT][MSTART] > 0 and best[ALIGNMENT][MEND] == len(best[ALIGNMENT][SEQ]):
        # -1 because originally had no +1 (whereas the above and below if statements had +1)
        delPos = best[ALIGNMENT][PTN].index('-') - 1
    # if a gap at the end ==>  deletion in the middle
    elif '-' in best[ALIGNMENT][PTN] and best[ALIGNMENT][PTN].index('-') + 1 < delPos:  # and best[1][4] < len(best[1][0]):
        delPos = best[ALIGNMENT][PTN].index('-')

    # find the location of mismatch
    misPos = -1
    # if it is Mismatched ==> length of alignment == length of pattern
    if len(best[ALIGNMENT][SEQ]) == len(patterns[bestInd][1]):
        misPos = 0
        while misPos < len(best[ALIGNMENT][SEQ]):
            # 5 is max score in the substitution matrix
            if subMatIUPAC[(best[ALIGNMENT][SEQ][misPos], patterns[bestInd][ALIGNMENT][misPos])] != 5:
                break
            misPos += 1
    # TODO: revise algorithm

    # 1-based
    """
    GGCCATCGGTCTCCCCC 
    [('alice', ('GGCCATCGGTCTCCCCC', 'GGTCACYG-TCTCYTCA', 43.0, 0, 16)),
     ('bob', ('--GG-CCATC-GGT-CTCCCCC', 'CAGGTBCAGCTGGTGCA-----', 31.0, 2, 16)),
     ('con', ('---GGCCATC-GGT-CTCCCCC', 'CARATGCAGCTGGTGCA-----', 21.0, 6, 16)),
     ('den', ('--GG-CCATC-GGT-CTCCCCC', 'SAGGTCCAGCTGGTACA-----', 31.0, 2, 16)),
     ('fur', ('GGCCATCGGTCTCCCCC-----', '---CA--GRTCACCTTGAAGGA', 26.0, 3, 14))]
    """
    if extend5end:
        return best[ID], misPos + 1, delPos + 1, extend5align(best[ALIGNMENT]), best[ALIGNMENT][MEND]
    # don't need to extend 5'end
    return best[ID], misPos + 1, delPos + 1, best[ALIGNMENT][MSTART], best[ALIGNMENT][MEND]


def extend5align(localAlignment):
    start = 0
    while start < len(localAlignment[1]) and localAlignment[1][start] == '-':
        start += 1
    return start if start > 0 else localAlignment[-2]


def splitFastaFile(fastaFile, totalFiles, seqsPerFile, filesDir,
                   prefix="", ext=".fasta", stream=None):
    if not exists(os.path.join(filesDir, prefix + "part" + str(int(totalFiles)) + ext)) and \
            not exists(os.path.join(filesDir, prefix + "part" + str(int(totalFiles)) + ".out")):
        # Split the FASTA file into multiple chunks
        printto(stream, "\tThe clones are distributed into multiple workers .. ")
        if not os.path.isdir(filesDir):
            os.makedirs(filesDir)
        if hasLargeMem():
            with safeOpen(fastaFile) as fp:
                recordsAll = SeqIO.to_dict(SeqIO.parse(fp, 'fasta'))
            queryIds = recordsAll.keys()
        else:
            # SeqIO.index can only open string filenames and they must be unzipped
            recordsAll = SeqIO.index(gunzip(fastaFile), 'fasta')
            # recordsAll.keys() is of type <dictionary-keyiterator object>, need to cast to list
            queryIds = list(recordsAll.keys())
        for i in range(totalFiles):
            ids = queryIds[i * seqsPerFile: (i + 1) * seqsPerFile]
            records = map(lambda x: recordsAll[x], ids)
            out = os.path.join(filesDir, prefix + 'part' + str(i + 1) + ext)
            SeqIO.write(records, out, 'fasta')


def writeSummary(filename, key, value):
    if os.path.exists(filename):
        with open(filename) as fp:
            string = fp.read()
    else:
        string = ""

    with open(filename, 'w') as fp:
        if re.search("^" + key + ":.*$", string, re.MULTILINE):
            string = re.sub("^" + key + ":.*$", "{}:{}".format(key, value), string, flags=re.MULTILINE)
        else:
            string += "{}:{}\n".format(key, value)
        fp.write(string)


def countSeqs(filename):
    # don't use line.startswith(">|@") here because
    # the quality score of FASTQ can start with @ too!
    _, ext = os.path.splitext(os.path.normpath(filename).replace(".gz", ""))
    ext = ext.lstrip(".")
    if ext in {"fasta", "fa"}:
        with safeOpen(filename) as fp:
            return sum([1 for _ in SeqIO.parse(fp, "fasta")])
    elif ext in {"fastq", "fq"}:
        with safeOpen(filename) as fp:
            return sum([1 for _ in SeqIO.parse(fp, "fastq")])
    else:
        raise ValueError("Unrecognized format {}, expected FASTA or FASTQ".format(ext))


def createIfNot(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def buildIgBLASTCommand(igdata, db, chain, species, domainSystem, blastInput, blastOutput, threads,
                        vOnly=False, protein=False, stream=None):
    if chain not in ['hv', 'lv', 'kv', 'klv']:
        raise ValueError("Unsupported chain type {}, expected one of 'lv', 'hv', 'kv', 'klv'.".format(chain))

    igdata = os.path.expandvars(igdata)
    db = os.path.expandvars(db)
    exe = IGBLASTP if protein else IGBLASTN
    # construct arguments for germline database (VDJ)
    germs = "V" if vOnly else "VDJ"

    cmd = {}
    for germ in germs:
        # if germline is D, then even the light chain will borrow imgt_<species>_ighd's database
        # otherwise it'll just be imgt_<species>_ig[kl(kl)][vj]
        chainLetter = 'h' if germ == 'D' else chain[:chain.find('v')]
        cmd['germline_db_' + germ] = os.path.join(db, 'imgt_{}_ig{}{}{}'.format(species,
                                                                                chainLetter,
                                                                                germ.lower(),
                                                                                '_p' if protein else ''))
    blast = ShortOpts(exe,
                      show_translation="",
                      extend_align5end="",
                      domain_system=domainSystem,
                      query=quote(blastInput),
                      organism=species,
                      auxiliary_data=os.path.join(igdata, 'optional_file', '{}_gl.aux'.format(species)),
                      outfmt=7,
                      num_threads=threads,
                      out=quote(blastOutput), **cmd)
    # printto(stream, "Executing : " + str(blast))
    return blast


def ntIUPACEqual(sequence, iupac):
    """
    compares 2 nucleotide sequence based on the IUPAC format. If they are not of equal length,
    they are automatically not equal. Comparison is case-insensitive

    ref: https://www.ddbj.nig.ac.jp/ddbj/code-e.html

    :param sequence: string
                normal ACGT nucleotide sequence
    :param iupac: string
                IUPAC format nucleotide sequence
    :return: bool
    """
    # string1 = unextended, # string2 = extended
    exMap = {
        'a': {'a'},
        'c': {'c'},
        'g': {'g'},
        't': {'t'},
        'm': {'a', 'c'},
        'r': {'a', 'g'},
        'w': {'a', 't'},
        's': {'c', 'g'},
        'y': {'c', 't'},
        'k': {'g', 't'},
        'v': {'a', 'c', 'g'},
        'h': {'a', 'c', 't'},
        'd': {'a', 'g', 't'},
        'b': {'c', 'g', 't'},
        'n': {'a', 'c', 'g', 't'}
    }
    if len(sequence) != len(iupac):
        return False
    for c1, c2 in zip(sequence, iupac):
        possibleMatches = exMap.get(c2.lower(), None)
        if None:
            raise ValueError("{} is not present under the IUPAC NT comparison mapping".format(c2))
        if c1.lower() not in possibleMatches:
            return False
    return True


def aaIUPACEqual(sequence, iupac):
    """
    compares 2 amino acid sequence based on the IUPAC format. If they are not of equal length,
    they are automatically not equal. Comparison is case-insensitive

    ref: https://www.ddbj.nig.ac.jp/ddbj/code-e.html

    :param sequence: string
                normal amino acid sequence
    :param iupac: string
                IUPAC amino acid sequence
    :return: bool
    """
    exMap = {
        'B': {'D', 'N'},
        'Z': {'Q', 'E'},
        'X': {'A', 'R', 'N', 'D', 'C', 'Q', 'E',
              'G', 'H', 'I', 'L', 'K', 'M', 'F',
              'P', 'O', 'S', 'U', 'T', 'W', 'Y',
              'V', 'B', 'Z', 'J'},
        'J': {'L', 'I'}
    }
    if len(sequence) != len(iupac):
        return False
    for c1, c2 in zip(sequence, iupac):
        if c1 != c2:
            possibleMatches = exMap.get(c2.upper(), None)
            if None:
                raise ValueError("{} is not present under the IUPAC AA comparison mapping".format(c2))
            if c1.upper() not in possibleMatches:
                return False
    return True

