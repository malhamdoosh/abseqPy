import os

SAMTOOLS_PROGRAM= 'samtools'
BGZIP = 'bgzip'
TABIX = 'tabix'
CLUSTALW = 'clustalw2'
CLUSTALOMEGA = 'clustal-omega'
WEBLOGO = 'weblogo'

tmp = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
MEM_GB = tmp/(1024.**3) 