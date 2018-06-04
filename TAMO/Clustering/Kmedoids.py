import string,sys,glob,random,pickle


'''
	Kmediods..py
	Currently contains a kMeans method
	carry out kMeans cluster using a precomputed distance matrix.
	Instead of using the mean of each cluster as the center, kMedoid uses 
	the central data point in the cluster

	kMeans clustering is non-deterministic, so it is advisable to use metaKMedoids_cluster()
	metaKMedoids_cluster() runs kMedoid clustering multiple times (determined by num_cycles)
	using different seeds (that can be supplied by the user)

	Both metaKMedoids_cluster() and kMedoids_cluster() return 3 items:
		1.  k indices representing the medoids
		2.  k-lists of indices representing the clusters, 
		3. the average distance to the medoid in each cluster
	

#Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
#All Rights Reserved
#
#Authors: Ernest Fraenkel & David Benjamin Gordon




	
'''

########################################################################
# "Meta-meta" Clustering Routines  ("Meta" and "Core" Routines below)  #
########################################################################

def JumpWithinKMedoids_cluster(distanceMatrix,kmax=10,num_cycles=5,max_kMedoids_iterations=1000,
			       min_dist=0.15,
			       seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	(There are better choices)  Increases k over a range (1-kmax) and selects k for which
	the change from k-1 "within_distance" is largest.
	'''
	kbest   = 0
	Wdists  = []
	Results = []
	for k in range(1,kmax,1):
		Wdist, medoids,clusters,distances = MinWithinKMedoids_cluster(
			distanceMatrix, k, num_cycles,
			max_kMedoids_iterations, seeds, verbose,
			print_fnctn, data,return_WD=1)
		Wdists.append(Wdist)
		Results.append((medoids,clusters,distances))
	diffs = []
	for i in range(len(Wdists)-1):
		delta = (Wdists[i] - Wdists[i+1])
		diffs.append(delta)
	print ['%5.1f'%x for x in Wdists]
	print '     ', ['%4.1f'%x for x in diffs]
	maxdelta = max(diffs)
	kbest    = diffs.index(maxdelta)+1
	print kbest
	print Results[kbest]
	sys.stdout.flush()
	return Results[kbest][0],Results[kbest][1],Results[kbest][2]


def bestKMedoids_cluster(distanceMatrix,kmax=10,num_cycles=5,max_kMedoids_iterations=1000,
			 min_dist=0.15,
			 seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	Descending k:  Stop Descending when distances between medoids becomes
	larger than min_dist
	'''
	
	kbest = kmax
	for k in range(kmax,0,-1):
		medoids,clusters,distances = metaKMedoids_cluster(
			distanceMatrix, k, num_cycles,
			max_kMedoids_iterations, seeds, verbose,
			print_fnctn, data)
		numzerodist = 0
		medidxs = medoids.values()
		dist = []
		for I in medidxs:
			for J in medidxs:
				if I >= J: continue
				dist.append(distanceMatrix[I][J])
				if distanceMatrix[I][J] < min_dist:
					numzerodist = numzerodist + 1
		if numzerodist == 0:
			kbest = k
			break
	#print "# Dist (unordered) [%s]"%(', '.join(['%6.4f'%x for x in dist]))
	try: print "# min,max dist:  %6.4f  %6.4f"%(min(dist),max(dist))
	except: pass
	return medoids,clusters,distances

def bestaveKMedoids_cluster(distanceMatrix,kmax=10,num_cycles=5,max_kMedoids_iterations=1000,
			    min_dist=0.15,
			    seeds=None,verbose=0, print_fnctn=None,data=None,kmin=0):
	'''
	[Ben\'s Current pick 12-02-03]
	Descending: Compute ave dist of members of a cluster to each of the other medoids.
	If any average distances are too small (<min_dist), reduce k
	''' 

	kbest = kmax
	for k in range(kmax,kmin,-1):
		#medoids,clusters,distances = metaKMedoids_cluster(
		medoids,clusters,distances =  MinWithinKMedoids_cluster(
			distanceMatrix, k, num_cycles,
			max_kMedoids_iterations, seeds, verbose,
			print_fnctn, data)
		numzerodist = 0
		avedists = []
		for medclust,medidx in medoids.items():
			for clust in range(k):              #Loop over clusters
				if clust == medclust: continue  #Skip my own cluster
				dtot = 0
				count = 0
				for J in clusters[clust]:
					#print clust,medidx,J,distanceMatrix[medidx][J]
					dtot = dtot + distanceMatrix[medidx][J]
					count = count + 1
				avedists.append(dtot/count)
		avedists.sort()
		if k!=1: print  "# K = %2d  mindist = %5.3f"%(k,avedists[0])
		sys.stdout.flush()
		#print "# Dist (unordered) [%s]"%(', '.join(['%5.3f'%x for x in avedists]))
		if k==1 or avedists[0] > min_dist:
			kbest = k
			break
	#print "# Dist (unordered) [%s]"%(', '.join(['%6.4f'%x for x in dist]))
	try: print "# min,max dist:  %6.4f  %6.4f"%(min(avedists),max(avedists))
	except: pass
	return medoids,clusters,distances

def metaminDiamKMedoids_cluster(distanceMatrix,kstart=1,num_cycles=5,max_kMedoids_iterations=1000,
			    min_dist=0.15,
			    seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	Ascending: increase k until the Diam (a.k.a. "distances") for each cluster are
	smaller than min_dist
	''' 

	kbest = kstart
	for k in range(kstart,len(distanceMatrix.keys())):
#		medoids,clusters,distances = metaKMedoids_cluster(
		medoids,clusters,distances =  minDiamKMedoids_cluster(
			distanceMatrix, k, num_cycles,
			max_kMedoids_iterations, seeds, verbose,
			print_fnctn, data)
		numzerodist = 0
		maxDiam = max(distances.values())
		print  "# K = %2d  maxDiam = %5.3f"%(k,maxDiam)
		sys.stdout.flush()
		#print "# Dist (unordered) [%s]"%(', '.join(['%5.3f'%x for x in avedists]))
		if maxDiam < min_dist:
			kbest = k
			break
	return medoids,clusters,distances

def bestupKMedoids_cluster(distanceMatrix,kmax=10,num_cycles=5,max_kMedoids_iterations=1000,
			 min_dist=0.15,
			 seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	Ascending (there are better choices)
	Ascending version of "best", which looks for the point at which the
	distances between medoids becomes too small
	'''
	kbest = kmax
	for k in range(2,kmax):
		medoids,clusters,distances = metaKMedoids_cluster(
			distanceMatrix, k, num_cycles,
			max_kMedoids_iterations, seeds, verbose,
			print_fnctn, data)
		numzerodist = 0
		medidxs = medoids.values()
		dist = []
		for I in medidxs:
			for J in medidxs:
				if I >= J: continue
				dist.append(distanceMatrix[I][J])
				if distanceMatrix[I][J] < min_dist:
					numzerodist = numzerodist + 1
		print "# Dist (unordered) [%s]"%(', '.join(['%6.4f'%x for x in dist]))
		print "Numzero: ",numzerodist,min(dist)
		if numzerodist > 0:
			kbest = k-1
			break
	medoids,clusters,distances = metaKMedoids_cluster(
		distanceMatrix, kbest, num_cycles,
		max_kMedoids_iterations, seeds, verbose,
		print_fnctn, data)
	numzerodist = 0
	medidxs = medoids.values()
	dist = []
	for I in medidxs:
		for J in medidxs:
			if I >= J: continue
			dist.append(distanceMatrix[I][J])
			if distanceMatrix[I][J] < min_dist:
				numzerodist = numzerodist + 1
	print "# Dist (unordered) [%s]"%(', '.join(['%6.4f'%x for x in dist]))
	print "Numzero: ",numzerodist
	#print "# Dist (unordered) [%s]"%(', '.join(['%6.4f'%x for x in dist]))
	try: print "# min,max dist:  %6.4f  %6.4f"%(min(dist),max(dist))
	except: pass
	return medoids,clusters,distances



#########################################################
# "Meta" Clustering Routines   ("Core" Routines below)  #
#########################################################


def MinWithinKMedoids_cluster(distanceMatrix,k,num_cycles=5,max_kMedoids_iterations=1000,seeds=None,verbose=0, print_fnctn=None,data=None,return_WD=None):
	'''
	Like original KMedoids_cluster, but the optimal iteration is picked based on 
	having the minimum 'within distance', which is computed by the withinDist
	function.
	WD = SUM_clusters (SUM_motifs (SUM_inter_motif_distances)))

	There is an optional flag to return WD to the parent, in case it wants to use it
	decide whether to repeat with a different value of k, as in "JumpWithin"

	'''
	if seeds and len(seeds)!= num_cycles:
		print 'num seeds != num cycles'
		sys.exit()
	if not seeds:
		seeds=[]
		for i in range(num_cycles):
			seeds.append(random.random())

	results=[]
	average_distances={}
	max_distances    ={}
	for i in range(num_cycles):
		#if verbose:	print 'run #',i
		medoids,clusters,distances=kMedoids_cluster(distanceMatrix,k,max_kMedoids_iterations,seeds[i])
		if verbose and print_fnctn :print_fnctn(clusters,data,distances,k)
		if not clusters: continue
		WD = withinDist(clusters,distanceMatrix)
		results.append((WD,medoids,clusters,distances))

	results.sort()
	WD,medoids,clusters,distances = results[0]
	if return_WD:
		return WD,medoids,clusters,distances
	else:
		return medoids,clusters,distances

def withinDist(clusters,data):
	#W = SUM_clusters (SUM_motifs (SUM_inter_motif_distances)))
	dtot = 0.0
	for clust,idxs in clusters.items():
		for i in idxs:
			for j in idxs:
				if i < j:
					dtot = dtot + data[i][j]
	return dtot
				


def metaKMedoids_cluster(distanceMatrix,k,num_cycles=5,max_kMedoids_iterations=1000,seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	###NOTE CHANGED FROM Fraenkel.Clustering###

	run kMedoid clustering num_cycles times
	return a three-tuple:
		1.  k indices representing the medoids
		2.  k-lists of indices representing the clusters, 
		3. the average distance to the medoid in each cluster
	if verbose =1 and a print function, print_fnctn, and the original data (not just the distance matrix) are supplied
	the program will print out intermediate clustering results
	'''
	if seeds and len(seeds)!= num_cycles:
		print 'num seeds != num cycles'
		sys.exit()
	if not seeds:
		seeds=[]
		for i in range(num_cycles):
			seeds.append(random.random())

	results=[]
	average_distances={}
	max_distances    ={}
	for i in range(num_cycles):
		if verbose:	print 'run #',i
		medoids,clusters,distances=kMedoids_cluster(distanceMatrix,k,max_kMedoids_iterations,seeds[i])
		if verbose and print_fnctn :print_fnctn(clusters,data,distances,k)
		results.append((medoids,clusters,distances))
		average_distances[i]=averageList(distances.values())
		max_distances[i]=max(distances.values())

	#items=average_distances.items()
	items=max_distances.items()
	items.sort(lambda x,y : cmp(x[1],y[1]))
	if verbose:
		print 'results of',num_cycles,'runs:'
		print items
		print 'best'
		print items[0][0]
	return results[items[0][0]]

def minDiamKMedoids_cluster(distanceMatrix,k,num_cycles=5,max_kMedoids_iterations=1000,
			    seeds=None,verbose=0, print_fnctn=None,data=None):
	'''
	'''
	if seeds and len(seeds)!= num_cycles:
		print 'num seeds != num cycles'
		sys.exit()
	if not seeds:
		seeds=[]
		for i in range(num_cycles):
			seeds.append(random.random())

	results=[]
	average_distances={}
	max_distances    ={}
	sumDiams         = {}
	for i in range(num_cycles):
		if verbose:	print 'run #',i
		medoids,clusters,distances=kMedoids_cluster(distanceMatrix,k,max_kMedoids_iterations,seeds[i])
		if verbose and print_fnctn :print_fnctn(clusters,data,distances,k)
		sum = 0
		for d in distances.values(): sum=sum+d
		results.append((medoids,clusters,distances))
		sumDiams[i] = sum
		average_distances[i]=averageList(distances.values())
		max_distances[i]=max(distances.values())

	#items=average_distances.items()
	items=sumDiams.items()
	items.sort(lambda x,y : cmp(x[1],y[1]))
	if verbose:
		print 'results of',num_cycles,'runs:'
		print items
		print 'best'
		print items[0][0]
	return results[items[0][0]]
	
			
	
##############################
# "Core" Clustering Routines #
##############################

def kMedoids_cluster(distanceMatrix,k,max_iterations=1000,seed ='12345'):
	'''
	use kMedoid clustering
	return a three-tuple:
		1.  k indices representing the medoids
		2.  k-lists of indices representing the clusters, 
		3. the average distance to the medoid in each cluster
	'''

	random.seed(seed)
	indices = distanceMatrix.keys()

	medoids={} #keys are =range(k); values are indices of distanceMatrix
	cluster_assignment={} #keys are the indices of distanceMatrix; values are selected from range(k)

	medoids=_initialize_medoids(indices,k)
	
	for iter in range(max_iterations):
		#print iter, 'iteration'
		#print 'medoids',medoids
	
		#old_clusters=cluster_assignment
		old_medoids=medoids

		#assign points to the closest medoid
		for cluster,index in medoids.items():
				cluster_assignment[index] = cluster
				#print index, 'is a medoid for', cluster
	
		for i in indices:
			if i in medoids.values():
				continue
			cluster_assignment[i] = _find_closest_medoid(distanceMatrix,medoids,i)
			#print i, 'assigned to', cluster_assignment[i]

		#calculate new medoids
			
		medoids = _calculate_medoids(distanceMatrix,cluster_assignment,k)
		#print 'medoids',medoids


		if old_medoids==medoids:
			break
	
	else:
		return None,None,None

	cluster={}
	avg_dist_in_cluster={}
	for i in range(k):
		cluster[i]=[]
	for i in indices:
		cluster[cluster_assignment[i]].append(i)
	for i in range(k):
		medoid=medoids[i]
		avg_dist_in_cluster[i]=0.
		for j in cluster[i]:
			avg_dist_in_cluster[i]+= distanceMatrix[j][medoid]
		avg_dist_in_cluster[i]/= len(cluster[i])

    	return medoids, cluster,avg_dist_in_cluster



def _initialize_medoids(indices,k,randomize=1):
	medoids={}
	if randomize:
    		random.shuffle(indices)
	for i in range(k):
		medoids[i]=indices[i]
	
    	return medoids
def  _find_closest_medoid(distanceMatrix,medoids,index):
	

	#print '_find_closest_medoid'
		
	closest_index = 0
    	closest_dist = distanceMatrix[index][medoids[0]]
    	for j in range(1,len(medoids)):
        	dist = distanceMatrix[index][medoids[j]]
        	if dist < closest_dist:
            		closest_dist = dist
            		closest_index = j
    	return closest_index

def _calculate_medoids(distanceMatrix,cluster_assignment,k):
	#print 'calculating new medoids'
	medoids={}
	clusters={}
	for i in range(k):
		clusters[i]=[]
	for i in cluster_assignment.keys():
		cluster_index=cluster_assignment[i]
		clusters[cluster_index].append(i)

	for cluster_index in clusters.keys():
		#print 'cluster', cluster_index, #clusters[cluster_index]
		cluster=clusters[cluster_index]
		cluster.sort()
		#print cluster
		sum_of_dist_dict ={}
		for i in cluster:
			sum_of_dist=0.
			for j in cluster:
				sum_of_dist+=distanceMatrix[i][j]
			sum_of_dist_dict[i]=sum_of_dist

		items=sum_of_dist_dict.items()
		items.sort(lambda x,y : cmp(x[1],y[1]))

		
		
		medoids[cluster_index]=items[0][0]
	return medoids
		






def print_list_clustering_results(clusters,data,avg_dist,k):
	print 'k=', k, 'average_distance=', averageList(avg_dist.values())
	for i in range(k):
		print
		print 'cluster', i
		print clusters
		for index in clusters[i]:
			print '\t',data[index]

def print_motif_clustering_results(clusters,data,avg_dist,k):
	print 'k=', k, 'average_distance=', averageList(avg_dist.values())
	for i in range(k):
		print
		print 'cluster', i
		for index in clusters[i]:
			motif= data[index]
			print '\t %-20s  %30s'%(motif.oneletter, motif.source)




def generate_random_1d (x=0,mult=1.,n=100):
	data =[]
	for i in range(n):
		data.append((x+mult*random.random()))
	return data
	
def generate_random_2d (x=0,y=0,mult=1.,n=100):
	data =[]
	for i in range(n):
		data.append((x+mult*random.random(), y+mult*random.random()))
	return data

def print_Dmat(D):
	Is = D.keys()
	Is.sort
	for I in Is:
		E = D[I]
		Js = E.keys()
		Js.sort()
		for J in Js:
			print '%4.2f '%E[J],
		print


def averageList(lst):
	if len(lst):
		return reduce(lambda x,y:x+y, lst) / float(len(lst))
	else:
		return -1
			

##################################################
# Alternative clustering methods                 #


	
if __name__ == '__main__':
	#import random

	#Data = [1,2,3,6,7,8,11,12,13]
	Data=[]
	Data+=generate_random_1d(10.,1,5)
	Data+=generate_random_1d(20.,1,5)
	Data+=generate_random_1d(30.,1,5)
	random.shuffle(Data)
	k=3


	
	distanceMatrix={}
	for i in range(len(Data)):
		distanceMatrix[i]={}
		for j in range(len(Data)):
			distanceMatrix[i][j] = abs(Data[i]-Data[j])
	
			

	medoids, clusters, avg_distance = metaKMedoids_cluster(distanceMatrix,k, verbose=0, print_fnctn=print_list_clustering_results,data=Data)
	for i in range(k):
			print 'cluster',i, 'avg dist', avg_distance
			for j in clusters[i]:
				print '\t',Data[j]

	

