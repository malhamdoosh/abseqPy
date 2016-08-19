import matplotlib as mpl
from collections import Counter
import math
mpl.use('Agg') # Agg

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from os.path import exists
from Bio import SeqIO
from numpy import Inf, mean, isnan
import random


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
    plt.close()
    
    
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
        weights = [1] * len(sizes)
    elif type(seqFile) == type([]):        
        sizes = map(lambda x: int(x) if not isnan(x) else 0, seqFile)
        weights = [1] * len(sizes)
    elif type(seqFile) == type(Counter()):
        sizes = seqFile.keys()
        weights = map(lambda x: seqFile[x], sizes)
    if removeOutliers:
        sizes, weights = excludeOutliers(sizes, weights)
    bins = max(sizes) - min(sizes) 
    if bins > maxbins:
        bins = bins / 2
    if maxbins == -1:
        bins = 20
    if bins == 0:
        bins = 1
#     print seqFile[:10], bins
    fig, ax = plt.subplots(figsize=(8,5))
    histcals, edges = np.histogram(sizes, bins = bins, range=autoscale,
                                   weights = weights, density=normed)
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
    plt.close()
    return histcals

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
    plt.close()
        
        
        
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
    plt.close()
    
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
    plt.close()
    
    
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
    plt.close()


def plotStatsHeatmap(data, sampleName, xyCol, axlabels, filename):    
    if (exists(filename)):
        print("File found ... " + filename.split('/')[-1])
        return
    x = data[xyCol[0]].tolist()
    y = data[xyCol[1]].tolist()
    total = len(x)
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
    title = 'Alignment Quality of Sample '
    title += '\nTotal is {:,}'.format(int(total)) 
    ax.set_title(title + sampleName)
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
    plt.close()
    
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
    


def excludeOutliers(values, weights, m =4.0):
    values = np.array(values)
    weights = np.array(weights)
    avg, std = weightedAvgAndStd(values, weights)
    sel =  abs(values - avg) <= m * std
    return values[sel].tolist(), weights[sel].tolist()





def weightedAvgAndStd(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))




    