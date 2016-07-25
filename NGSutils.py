import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from os.path import exists
import os
import sys
from Bio import SeqIO, AlignIO
from pandas.core.frame import DataFrame
from numpy import Inf, mean, isnan
from Bio.SeqRecord import SeqRecord
from config import CLUSTALOMEGA
from Bio.Align.Applications._Clustalw import ClustalwCommandline
from Bio.Seq import Seq
from collections import Counter
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import random
import re

def fastq2fasta(fastqFile, outputDir):
    # FASTQ to FASTA
# awk 'NR % 4 == 1 {print ">" $0 } NR % 4 == 2 {print $0}' my.fastq > my.fasta
    filename = fastqFile.split('/')[-1]    
    seqOut = outputDir + "seq/"
    if (not os.path.isdir(seqOut)):
        os.system("mkdir " + seqOut)
    filename = seqOut + filename.replace(filename.split('.')[-1], 'fasta')
    if exists(filename):
        print ("\tThe FASTA file was found!")
        return filename
    print(fastqFile + " is being converted into FASTA ...")
    command = ("awk 'NR % 4 == 1 {sub(\"@\", \"\", $0) ; print \">\" $0} NR % 4 == 2 "
               "{print $0}' " + fastqFile + " > " + filename
               )
    
    os.system(command)
    return filename 

def runIgblastn(blastInput, chain, threads = 8, db='$IGBLASTDB'):
    # Run igblast on a fasta file  
      
    blastOutput = blastInput.replace('.' + blastInput.split('.')[-1], '.out')
    if (exists(blastOutput)):
        print("\tBlast results were found ... " + blastOutput.split("/")[-1])
        return blastOutput 
    print('\tRunning igblast ... ' + blastInput.split("/")[-1])
    if (chain == 'hv'):
        command = ("igblastn -germline_db_V " + db+"/imgt_human_ighv -germline_db_J " 
                   "" + db+"/imgt_human_ighj -germline_db_D " + db+"/imgt_human_ighd -domain_system imgt "
                   "-query %s -organism human -auxiliary_data optional_file/human_gl.aux "  
                   "-show_translation -outfmt 7 -num_threads %d -out %s"
                   )
    elif (chain == 'kv'):
        command = ("igblastn -germline_db_V " + db+"/imgt_human_igkv -germline_db_J " 
                   "" + db+"/imgt_human_igkj -germline_db_D " + db+"/imgt_human_ighd -domain_system imgt "
                   "-query %s -organism human -auxiliary_data optional_file/human_gl.aux "  
                   "-show_translation -outfmt 7 -num_threads %d -out %s"
                   )
    elif (chain == 'lv'):
        command = ("igblastn -germline_db_V " + db+"/imgt_human_iglv -germline_db_J " 
                   "" + db+"/imgt_human_iglj -germline_db_D " + db+"/imgt_human_ighd -domain_system imgt "
                   "-query %s -organism human -auxiliary_data optional_file/human_gl.aux "  
                   "-show_translation -outfmt 7 -num_threads %d -out %s"
                   )
    else:
        print('ERROR: unsupported chain type.')     
        sys.exit()   
        
    os.system(command % (blastInput, threads, blastOutput))
    return blastOutput

def runIgblastp(blastInput, chain, threads = 8, db='$IGBLASTDB'):
    # Run igblast on a fasta file        
    blastOutput = blastInput.replace('.' + blastInput.split('.')[-1], '.out')
    if (exists(blastOutput)):
        print("\tBlast results were found ... " + blastOutput.split("/")[-1])
        return blastOutput 
    print('\tRunning igblast ... ' + blastInput.split("/")[-1])
    if (chain == 'hv'):
        command = ("igblastp -germline_db_V " + db+"/imgt_human_ighv_p " 
                   "-domain_system imgt "
                   "-query %s -organism human "  
                   "-outfmt 7 -num_threads %d -out %s"
                   )
    elif (chain == 'kv'):
        command = ("igblastp -germline_db_V " + db+"/imgt_human_igkv_p " 
                   "-domain_system imgt "
                   "-query %s -organism human "  
                   "-outfmt 7 -num_threads %d -out %s"
                   )
    elif (chain == 'lv'):
        command = ("igblastp -germline_db_V " + db+"/imgt_human_iglv_p " 
                   "-domain_system imgt "
                   "-query %s -organism human "  
                   "-outfmt 7 -num_threads %d -out %s"
                   )
    else:
        print('ERROR: unsupported chain type.')     
        sys.exit() 
        
    
    os.system(command % (blastInput, threads, blastOutput))
    return blastOutput


def writeCountsToFile(dist, filename):
    # This function prints the distribution counts into a text file
    with open(filename, 'w') as out:
        out.write('IGHV Class, Count, Proportion \n')
        total = sum(dist.values()) * 1.0
        for k in sorted(dist, key=dist.get, reverse=True):
            out.write(k + ',' + `dist[k]` + ',' + ("%.2f" % (dist[k] / total * 100)) + '\n')
        out.write('TOTAL, ' + `total` + ', 100 ')
    print("A text file has been created ... " + filename)
            
def plotSeqLenDistClasses(seqFile, sampleName, outputFile, fileFormat='fasta', maxLen=Inf):
    if (exists(outputFile)):
        print("File found ... " + outputFile.split('/')[-1])
        return
    print("The sequence length distribution of each gene family is being calculated ...")
    ighvDist = {}
    ighvSizes = {}
    for rec in SeqIO.parse(seqFile, fileFormat):
        if (len(rec) <= maxLen):
            if (rec.id.split('|')>1):
                ighvID = rec.id.split('|')[1]
            else:
                ighvID = rec.id
            id = ighvID.split('-')[0].split('/')[0]
            if ighvDist.get(id, None) is None:
                ighvDist[id] = 0
                ighvSizes[id] = []
            ighvSizes[id].append(len(rec))
            ighvDist[id] += 1
 
    plotDist(ighvDist, sampleName, outputFile)
    # box plot of sequence length in each class
    fig, ax = plt.subplots()
    classes = sorted(ighvDist, key=ighvDist.get, reverse=True)
    ax.boxplot(map(lambda x:ighvSizes[x], classes))
    ind = np.arange(1, len(classes)+1)
    ax.set_xticks(ind)
    ax.set_xticklabels(classes, rotation=45)    
    ax.set_title("Sequence Lengths in " + sampleName)
    outputFile = '/'.join(outputFile.split('/')[:-1] + ["box_"+outputFile.split('/')[-1]])
    fig.savefig(outputFile, dpi=300)
    for k in classes:
        print(k, ighvDist[k], min(ighvSizes[k]), max(ighvSizes[k]))
    
def excludeOutliers(list, m =4.0):
    data = np.array(list)
    return data[abs(data - np.mean(data)) <= m * np.std(data)].tolist()
    
def plotSeqLenDist(seqFile, sampleName, outputFile, fileFormat='fasta', 
                   maxLen=Inf, histtype='bar', dna=True, 
                   autoscale=None, maxbins=20, seqName='', normed=False,
                   removeOutliers=False):
    if (exists(outputFile)):
        print("File found ... " + outputFile.split('/')[-1])
        return
    print("The sequence length distribution is being calculated ...")
    if (type("") == type(seqFile)):
        sizes = [len(rec) for rec in SeqIO.parse(seqFile, fileFormat) if len(rec) <= maxLen]
    else:
        assert type(seqFile) == type([])
        sizes = map(lambda x: int(x) if not isnan(x) else 0, seqFile)
    if removeOutliers:
        sizes = excludeOutliers(sizes)
    bins = max(sizes) - min(sizes) 
    if bins > maxbins:
        bins = bins / 2
    if maxbins == -1:
        bins = 20
    if bins == 0:
        bins = 1
#     print seqFile[:10], bins
    fig, ax = plt.subplots(figsize=(8,5))
    histcals, edges = np.histogram(sizes, bins = bins, range=autoscale, density=normed)
    binWidth = edges[1] - edges[0] 
    ax.bar(edges[:-1], histcals * binWidth, binWidth)
#     histcals = ax.hist(sizes, bins=bins, histtype=histtype, range=autoscale,
#                        normed=normed)
    title = "{:,} Sequences {} in {} \nLengths {:d} to {:d}"
    ax.set_title(title.format(len(sizes), 'of ' + seqName if seqName!='' else '',
                               sampleName, min(sizes), max(sizes)))
    if dna:
        ax.set_xlabel("Sequence Length (bp)")
    else:
        ax.set_xlabel("Sequence Length (aa)")
    if autoscale:
        ax.set_xticks(np.arange(autoscale[0], autoscale[1]+1, 5))
    if (not normed):
        ax.set_ylabel("Count")
#         ax.set_ylim(top=len(sizes))
    else:
        ax.set_ylabel("Proportion")
#         ax.set_ylim(top=1)
    fig.savefig(outputFile, dpi=300)
    return histcals

def findBestAlignment(seq, query, dna=False, offset=0, show=False):
    if not dna:
        alignments = align.localds(seq.replace('*', 'X'), query, matlist.blosum62,-100, -100)
    else:
        alignments = align.localms(seq, query, 1,-2, -2, -2)
    
#     print(seq, query, alignments)
    scores = [a[2] for a in alignments]
    if (len(scores) == 0):
#         print(seq, query, alignments)
#         raise
        return -1, -1
    best = scores.index(max(scores))
    if show:
        print(format_alignment(*alignments[best]))
    # return alignment start and end
    return int(offset + alignments[best][-2] + 1), int(offset + alignments[best][-1]) # 1-based

'''
    Extract a protein fragment from a protein sequence based on DNA positions
'''
def extractProteinFrag(protein, start, end, offset=0, trimAtStop=False):
    if (np.isnan(start) or np.isnan(end)):
        return None
    # start and end are 0-based positions
    try:
        if (start != -1):
            s = int(round((start - offset - (start-offset)%3) / 3.))
        else:
            s = 0
        if end != -1:
            e = int(round((end - offset- (start-offset)%3) / 3.))
        else:
            e = len(protein)
        if s < e:        
            frag = protein[s:e]
        elif s == e:
            frag = protein[e-1]
        else:
            return None
        if (frag is None):
            print("Extract Protein Fragment",protein, start, end)
            return None
        if trimAtStop and ('*' in frag):
            frag = frag[:frag.index('*')]
        return frag
    except:
        print("Extract Protein Fragment",protein, start, end)
        return None
    
def plotSeqDuplication(seqs, filename, labels, title='', grouped=False):
    if (exists(filename)):
        print('File found ... ' + filename.split('/')[-1])
        return
    print("The duplication of VH sequences is being estimated .... ")
    fig, ax = plt.subplots(figsize=(8,5))
    ax.grid()
    ax.set_xlabel('Duplication Level')
    ax.set_ylabel('Proportion of Duplicated Sequences')
    if (not grouped):
        ax.set_title(title + '\nTotal is {:,}'.format(int(len(seqs[0]))))
    else:
        ax.set_title(title + '\nTotal is {:,}'.format(sum(map(lambda x: len(x), seqs))))
    for setSeqs, l in zip(seqs, labels):
        dup = {}
        for s in setSeqs:
            dup[s] = dup.get(s, 0) + 1
#             print(dup[s])
#         print("here1")
        freqs = dup.values()
        freqs.sort()
#         print("here2", len(freqs), freqs[0], freqs[-1])
        ticks = np.linspace(10, 10000, 100).tolist()
#         ticks.append(freqs[-1])
        y = []
        freqs = np.array(freqs)
        for x in ticks:
            y.append(sum(1.0 * freqs[freqs >= x]) / len(setSeqs))
        # scale ticks to [10, 20]
        ticks = map(lambda x : (x - 10) * (20 - 10) / (10000 - 10) + 10 , ticks)
        less10Ticks = []
        less10Y = []
        for i in range(1, 10):
            less10Ticks.append(i)
            less10Y.append(sum(freqs == i) * 1.0 / len(setSeqs))        
        y = less10Y + y
        ticks = less10Ticks + ticks 
#         print("here3", len(y), ticks[:20], y[:20])
        ax.plot(ticks,y , label = l)
#         break
    xticks = range(1, 21, 2)
    ax.set_xticks(xticks)
    xlabels = range(1, 10, 2) 
    xlabels += map(lambda x: '>' + (`int(x) - int(x)%100` if x > 100 else `int(x)`), 
                   np.linspace(10, 10000, (len(xticks) - len(xlabels)) * 2 ).tolist()[1::2])
    ax.set_xticklabels(xlabels)
    ax.legend()
    fig.savefig(filename, dpi=300)
        

def plotSeqDiversity(seqs, filename, labels, title=''):
    if (exists(filename)):
        print('File found ... ' + filename.split('/')[-1])
        return
    print("The diversity of VH sequences is being estimated .... ")
    fig, ax = plt.subplots(figsize=(8,5))
    ax.grid()
    ax.set_xlabel('No of Sequences')
    ax.set_ylabel('Number of Deduplicated Sequences')
    ax.set_title(title)

    for setSeqs, l in zip(seqs, labels):
        ticks = []
        S = 10
        while S < len(setSeqs):
            ticks.append(S)
            S = int(S * 1.5)
        ticks.append(len(setSeqs))
#         print ticks
        pt = []
        for j in ticks:
            hs = [ len(set(random.sample(setSeqs, j))) for k in range(5) ]
            pt.append((j, hs))
        pt.append((len(setSeqs), len(set(setSeqs))))
        ax.plot([ d[0] for d in pt ], [ mean(d[1])*1.0 for d in pt ], label = l)
    ax.legend()
    fig.savefig(filename, dpi=300)
    
'''
    Plot Venn diagrams using the matplotlib_venn package
'''
def plotVenn(sets, filename):
    if (exists(filename)):
        print("File found ... " + filename.split('/')[-1])
        return
    fig, ax = plt.subplots()
    if (len(sets) == 2):
        from matplotlib_venn  import venn2
        venn2(sets.values(), sets.keys())
    elif (len(sets) == 3):
        from matplotlib_venn  import venn3
        venn3(sets.values(), sets.keys())
    else:
        raise
    fig.savefig(filename, dpi=300)
    
    
def plotDist(ighvDistfam, sampleName, filename, title='', proportion=True, 
             rotateLabels=True, vertical=True, sortValues=True, top=15):   
    if (exists(filename)):
        print("File found ... " + filename.split('/')[-1])
        return
    # This function creates bar plot for the distribution counts/proportions 
    if sortValues:
        classes = sorted(ighvDistfam, key=ighvDistfam.get, reverse=True)
    else:
        classes = ighvDistfam.keys()
        classes.sort()
    if (len(classes) > top):
        classes = classes[:top]
    if not vertical:
        classes = classes[::-1]
    total = sum(ighvDistfam.values()) * 1.0
#     if (proportion):
    stats = map(lambda x: ighvDistfam[x] / total * 100, classes)
#     else:
    if (len(stats) < 1):
        return
    ind = np.arange(len(classes))
    
    fig, ax = plt.subplots(figsize=(8,5) if vertical else (5,8))
    ax.grid()
    if len(classes) > 10 :
        width = 0.4        
    else:
        width = 0.6
    if proportion:        
        topvalFormat = '{:.2f}'
    else:
        topvalFormat = '{:,}'
    # Create the bar plot and format it      
    if vertical:
        rects = ax.bar(ind, stats, width)
        ax.set_xticks(ind+width/2)
        ax.set_ylim(top=max(stats) * 1.1)
        if rotateLabels:
            ax.set_xticklabels(classes, rotation=45)
        else:
            ax.set_xticklabels(classes)
        ax.set_ylabel('Proportion (%)')
        # write the proportion on the top of each bar
        for rect in rects:
            height = rect.get_height()
            if not proportion:
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, 
                    (topvalFormat.format(int(np.round(height * total / 100)))),
                    ha='center', va='bottom', size=10, color='red')
            else:
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, 
                    (topvalFormat.format(height)),
                    ha='center', va='bottom', size=10, color='red')
    else:
        rects = ax.barh(ind, stats, width)
        ax.set_yticks(ind+width/2)
        ax.set_xlim(right=max(stats) * 1.1)
        if rotateLabels:
            ax.set_yticklabels(classes, rotation=45)
        else:
            ax.set_yticklabels(classes)
        ax.set_xlabel('Proportion (%)')
        # write the proportion on the top of each bar
        for rect in rects:
            width = rect.get_width()
            if not proportion:                
                ax.text(0.8 + width, rect.get_y()+rect.get_height()/2.,  
                    (topvalFormat.format(int(np.round(width * total / 100)))),
                    ha='center', va='bottom', size=10, color='red')
            else:
                ax.text(0.8 + width, rect.get_y()+rect.get_height()/2.,  
                    (topvalFormat.format(width)),
                    ha='center', va='bottom', size=10, color='red')
            
    
    if (title == ''):
        title = 'IGV Abundance in Sample ' + sampleName 
    title += '\nTotal is {:,}'.format(int(total)) 
    ax.set_title(title)  
    plt.tight_layout()
    
   
    
    fig.savefig(filename, dpi=300)


def plotStatsHeatmap(data, sampleName, xyCol, axlabels, filename):    
    x = data[xyCol[0]].tolist()
    y = data[xyCol[1]].tolist()
    BINS = 10
#     fig, ax = plt.subplots()        
#     ax.scatter(x, y, s=3, alpha=0.5, edgecolors='none' )
#     ax.set_xlabel(axlabels[0])   
#     ax.set_ylabel(axlabels[1])    
#     ax.set_title('Alignment Quality of Sample ' + sampleName)       
#     fig.savefig(filename,
#                 dpi=300)
    
    # plot as heatmap
    fig, ax = plt.subplots() 
    heatmap, xedges, yedges = np.histogram2d(x,y, bins=BINS)   
#     print xedges
#     print yedges
    heatmap = heatmap / np.sum(heatmap) * 100
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#     c = cmap_discretize('coolwarm', 5)
    cax = ax.imshow(heatmap.transpose()[::-1], cmap='jet', interpolation='nearest', 
                    extent=extent)
    ax.set_xlabel(axlabels[0])   
    ax.set_ylabel(axlabels[1])    
    ax.set_title('Alignment Quality of Sample ' + sampleName)
    ax.set_xticks(np.array(xedges).astype(int))
    ax.set_yticks(np.array(yedges).astype(int))
    ax.tick_params(axis='both', which='major', labelsize=8)
#     ax.set_xticklabels(np.round(np.linspace(xedges[0], xedges[-1], BINS/2)))
    cbar = fig.colorbar(cax,    
                        ticks = np.linspace(np.min(heatmap),
                                            np.max(heatmap),
                                            5),         
                orientation='horizontal')
#     print(np.percentile(heatmap, [0, 25, 50, 75, 100]))
#     labels = np.percentile(heatmap, [0, 25, 50, 75, 100])        
#     cbar.set_ticklabels(labels)

    forceAspect(ax,aspect=1)
    fig.savefig(filename, dpi=300)
    
def mergeReads(readFile1, readFile2, threads=3, merger='leehom', outDir="./"):    
    readFile = readFile1.split("/")[-1]
    outputPrefix = readFile.replace("_" + readFile.split('_')[-1], '')
    seqOut = outDir + "seq/"
    if (not os.path.isdir(seqOut)):
        os.system("mkdir " + seqOut)
    mergedFastq = seqOut 
    if (merger == 'pear'):        ### MERGE using PEAR            
        mergedFastq += outputPrefix + '.assembled.fastq'
        if (not exists(mergedFastq)):
            print("%s and %s are being merged ..." % (readFile1.split('/')[-1]
                                              , readFile2.split('/')[-1])) 
            command = "pear -f %s -r %s -o %s -j %d -v 15 -n 350"
            os.system(command % (readFile1, readFile2, outputPrefix, threads))
            os.system("mv %s.* %s" % (outputPrefix, seqOut))            
        else:
            print(".../" + mergedFastq.split("/")[-1] + ' was found!')
    elif (merger == 'leehom'):        
        mergedFastq += outputPrefix + '.fq'
        if (not exists(mergedFastq)):
            print("%s and %s are being merged ..." % (readFile1.split('/')[-1]
                                              , readFile2.split('/')[-1])) 
            command = "leeHom -fq1 %s -fq2 %s -fqo %s --ancientdna --verbose"
            os.system(command % (readFile1, readFile2, outputPrefix))
            os.system('gunzip ' + mergedFastq + '.gz')
            os.system("mv %s.* %s" % (outputPrefix, seqOut))
            os.system("mv %s_r* %s" % (outputPrefix, seqOut))
        else:
            print(".../" + mergedFastq.split("/")[-1] + ' was found!')
    elif (merger == 'flash'):        
        mergedFastq += outputPrefix + '.extendedFrags.fastq'
        if (not exists(mergedFastq)):
            print("%s and %s are being merged ..." % (readFile1.split('/')[-1]
                                              , readFile2.split('/')[-1])) 
            # the merger params souldn't be hardcoded
            command = "flash %s %s -t %d -o %s -r 300 -f 500 -s 150"            
            os.system(command % (readFile1, readFile2, threads, outputPrefix))
            os.system("mv %s.* %s" % (outputPrefix, seqOut))            
        else:
            print(".../" + mergedFastq.split("/")[-1] + ' was found!')
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
        raise Exception("Uknowne short reads merger is selected")
        
    return os.path.abspath(mergedFastq)


def writeTableIntoFile(table, filename):
    df = DataFrame(table)
    df.to_csv(filename, sep='\t', header=True, index=True)    
    print("Text file has been written to " + filename)

def writeListToFile(items, filename):
    out = open(filename, 'w')
    out.write("\n".join(items))
    out.close()

def loadIGVSeqsFromFasta(filename):
    ighvSeqs = {}
    for rec in SeqIO.parse(filename, 'fasta'): 
        ighv = rec.id.split('|')[1].strip()
        if (ighvSeqs.get(ighv, None) is None):
            ighvSeqs[ighv] = []
        ighvSeqs[ighv].append(str(rec.seq))
        
    return ighvSeqs
                
def compressSeqGeneLevel(seqDict):    
    geneLevel = {}
    for ighv in seqDict.keys():
        gene = ighv.split('*')[0]
        if (geneLevel.get(gene, None) is None):
            geneLevel[gene] = []
        geneLevel[gene] += seqDict[ighv]
    return geneLevel

def compressSeqFamilyLevel(seqDict):    
    familyLevel = {}
    for ighv in seqDict.keys():
        fam = ighv.split('-')[0].split('/')[0]
        if (familyLevel.get(fam, None) is None):
            familyLevel[fam] = []
        familyLevel[fam] += seqDict[ighv]
    return familyLevel

def compressCountsGeneLevel(countsDict):
    geneLevel = Counter()
    for k in countsDict.keys():
        ksub = k.split('*')[0]
        geneLevel[ksub] = geneLevel.get(ksub, 0) + countsDict[k]
    return geneLevel

def writeCountsCategoriesToFile(countsVariant, sampleName, filePrefix, title=''):
    writeCountsToFile(countsVariant,
            filePrefix + 'variant.csv')
    # gene level
    countsVariant = compressCountsGeneLevel(countsVariant)
    writeCountsToFile(countsVariant,
        filePrefix + 'gene.csv')
    plotDist(countsVariant, sampleName, 
             filePrefix + 'gene.png', 
             title)
    # family level
    countsVariant = compressCountsFamilyLevel(countsVariant)
    writeCountsToFile(countsVariant,
        filePrefix + 'family.csv')
    plotDist(countsVariant, sampleName, 
             filePrefix + 'family.png',
             title)
            
def compressCountsFamilyLevel(countsDict):
    familyLevel = Counter()
    for k in countsDict.keys():
        ksub = k.split('-')[0].split('/')[0].rstrip('D')
        familyLevel[ksub] = familyLevel.get(ksub, 0) + countsDict[k]
    return familyLevel
            
'''
    perform multiple sequence alignment using CLUSTAL
'''
def alignListOfSeqs(signals):
    L = map(len, signals)
    print("\t\t%d sequences are being aligned using CLUSTAL-OMEGA (L in [%d, %d])... " % (len(L), min(L), max(L)))
    tempSeq = "csl_temp_seq.fasta"
    tempAlign = tempSeq.replace('.fasta', '.aln')
    seqs = []
    for i in range(len(signals)):
        seqs.append(SeqRecord(Seq(signals[i]), id='seq' + `i`))
    SeqIO.write(seqs, tempSeq, 'fasta')
    clustalw = ClustalwCommandline(CLUSTALOMEGA, infile=tempSeq, 
                    outfile=tempAlign)
    stdout, stderr = clustalw() 
    
    alignment = AlignIO.read(tempAlign, 'clustal')
    alignedSeq = []
    for rec in alignment:
        alignedSeq.append(str(rec.seq))
    os.system("rm %s %s " % (tempSeq, tempAlign) )
    return alignedSeq

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
            iupacSubMat[(colHeads[i], row[0])] = float(row[i+1])
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
        
'''
    A function to find the best matched pattern in a list of patterns
    and classify the type of the alignment (intact, indelled, mismatched, unknown)
'''            
def findBestMatchedPattern(seq, patterns):
    scores = []
    # align the sequence against all possible patterns
    for (id, pattern, maxScore) in patterns:
#         print(seq, pattern)
        alignments = align.localds(seq.upper(), pattern, subMatIUPAC, -5, -5)
        if (len(alignments) > 1):
            localScores = [a[2] for a in alignments]
            alignment = alignments[localScores.index(max(localScores))]
        elif (len(alignments) > 0):
            alignment = alignments[0]
        else:
            raise Exception("Couldn't align ... ", seq, patterns)            
        if alignment:
            alignLen = alignment[-1] - alignment[-2]
            scores.append((id, alignment))
            ## if the sequence exactly matches one of the patterns ==> intact
            if (alignment[2] == maxScore and 
                alignLen == len(pattern) and 
                '-' not in alignment[0] and
                '-' not in alignment[1]):
                return (scores[-1][0], "Intact", 0)
        else:
            scores.append((id, ('', '', 0)))
    # if no exact matching ==> find the best alignment (pattern)
    if (len(scores) > 1):       
        tmp = map(lambda x : x[1][2], scores)
        bestInd = tmp.index(max(tmp))       
    elif len(scores) == 1:
        bestInd = 0
    else:
        return ("None", "Unknown", 0)
    # classify the allignment type ==> insertion, deletion, mismatches
    best = list(scores[bestInd])    
    best[1] = list(best[1])
    if best[1][2] == 0:
        return ("None", "Unknown", 0)
    # Find the position of Indel/Mismatch
    # remove starting indels
    if best[1][1].startswith('-'):
        i = 0
        while best[1][1][i] == '-':
            i += 1
        best[1][0] = best[1][0][i:]
        best[1][1] = best[1][1][i:]
    # find the location of insertion or deletion 
    delPos = -1
    if '-' in best[1][0]:
        delPos=  best[1][0].index('-') + 1
    # if there is a gap at the beginning ==> happened because of insertion/deletion in the middle
    if '-' in best[1][1] and  best[1][1].index('-') > delPos and best[1][3] > 0 and best[1][4] == len(best[1][0]):
        delPos = best[1][1].index('-')
    # if a gap at the end ==>  deletion in the middle
    elif '-' in best[1][1] and  best[1][1].index('-') + 1 < delPos: # and best[1][4] < len(best[1][0]):
        delPos = best[1][1].index('-') + 1
        
    if delPos != -1:    
        return (best[0], "Indelled", delPos) # 1-based
    else:
        # if it is Mismatched ==> length of alignment == length of pattern
        try:           
            assert len(best[1][0]) == len(patterns[bestInd][1])
            misPos = 0
            while (misPos < len(best[1][0])):
                # 5 is max score in the substitution matrix
                if (subMatIUPAC[(best[1][0][misPos], patterns[bestInd][1][misPos])] != 5):
                    break
                misPos += 1 
        except:
            raise Exception("Unexpected behaviour:" + seq +  " " + scores + " " + best)            
        return (best[0], "Mismatched", misPos+1) # 1-based
     
        
            
def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in xrange(N+1) ]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)