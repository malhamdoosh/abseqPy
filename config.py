'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''


import os
import sys

#SAMTOOLS_PROGRAM= 'samtools'
#BGZIP = 'bgzip'
#TABIX = 'tabix'
VERSION = '1.1.4'
CLUSTALW = 'clustalw2'
CLUSTALOMEGA = 'clustalo'
WEBLOGO = 'weblogo'
FASTQC = 'fastqc'
DEFAULT_MERGER = 'leehom'
DEFAULT_TOP_CLONE_VALUE = 100
RSCRIPT_PAIRING_SEPARATOR = ';'
RSCRIPT_SAMPLE_SEPARATOR = '|'

# consensus protein of HV http://discovery.ucl.ac.uk/15808/1/15808.pdf
FR4_CONSENSUS = {'hv':"WGQGTXVTVSS", 'kv':'FGGGTQ', 'lv':'FGGGTQ'}
FR4_CONSENSUS_DNA = {'hv':"TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC", 
                     'kv':'TTCGGCGGCGGCACCCAG', 'lv':'TTCGGCGGCGGCACCCAG'}

GB = (1024.**3)

# sorry darwin people, you need psutil because sysconf can't locate 'sc_phys_pages'
if sys.platform == 'darwin':
    from psutil import virtual_memory
    mem = virtual_memory()
    MEM_GB = mem.total/GB
else:
    tmp = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
    MEM_GB = tmp/GB
