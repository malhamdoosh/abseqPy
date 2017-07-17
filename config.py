'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import os

#SAMTOOLS_PROGRAM= 'samtools'
#BGZIP = 'bgzip'
#TABIX = 'tabix'
VERSION = '1.1.1'
CLUSTALW = 'clustalw2'
CLUSTALOMEGA = 'clustal-omega'
WEBLOGO = 'weblogo'
FASTQC = 'fastqc'

# consensus protein of HV http://discovery.ucl.ac.uk/15808/1/15808.pdf
FR4_CONSENSUS = {'hv':"WGQGTXVTVSS", 'kv':'FGGGTQ', 'lv':'FGGGTQ'}
FR4_CONSENSUS_DNA = {'hv':"TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC", 
                     'kv':'TTCGGCGGCGGCACCCAG', 'lv':'TTCGGCGGCGGCACCCAG'} 

tmp = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
MEM_GB = tmp/(1024.**3) 
