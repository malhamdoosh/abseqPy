'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''


import os
import sys

# ==========================================
#           ABSEQ's VERSION
# ==========================================
ABSEQROOT = os.path.abspath(os.path.dirname(__file__))
VERSION = '1.1.6'

# ====================================================================================
#           ABSEQ's EXTERNAL DEPENDENCIES
# ====================================================================================
# NOTE TO PROGRAMMER: IF YOU CHANGE 3rd_party TO SOME OTHER DIRECTORY NAME, MAKE SURE YOU CHANGE
# IT IN setup.py AND MANIFEST.in TOO! (just search for this comment and you'll find the exact location)
EXTERNAL_DEP_DIR = '3rd_party/'
# 1. clustal omega
CLUSTALOMEGA = 'clustalo'
# 2. FASTQC
FASTQC = 'fastqc'
# 3. leeHom
LEEHOM = 'leeHomMulti'
# 4.a IgBlastn
IGBLASTN = 'igblastn'
# 4.b IgBlastp
IGBLASTP = 'igblastp'


# temporarily overrides PATH variable with EXTERNAL_DEP_DIR/bin, IGBLASTDB and IGDATA (if they exist)
class PriorityPath:
    def __init__(self):
        self.updated = False
        self.old_env = os.environ.copy()
        _env = os.environ.copy()

        # if the BIN directory exists, append it to the front of PATH variable
        override_path = os.path.abspath((ABSEQROOT + '/' + EXTERNAL_DEP_DIR + '/bin/').replace('//', '/'))
        if os.path.exists(override_path):
            _env['PATH'] = override_path + os.pathsep + _env['PATH']
            self.updated = True

        # if the igdata dir exists, override it irrespective of if there's already a IGDATA env
        override_igdata = os.path.abspath((ABSEQROOT + '/' + EXTERNAL_DEP_DIR + '/igdata/').replace('//', '/'))
        if os.path.exists(override_igdata):
            _env['IGDATA'] = override_igdata
            self.updated = True

        # if the igdb dir exists, override it irrespective of if there's already a IGBLASTDB env
        override_igdb = os.path.abspath((ABSEQROOT + '/' + EXTERNAL_DEP_DIR + '/databases/').replace('//', '/'))
        if os.path.exists(override_igdb):
            _env["IGBLASTDB"] = override_igdb
            self.updated = True

        if self.updated:
            os.environ.clear()
            os.environ.update(_env)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.updated:
            os.environ.clear()
            os.environ.update(self.old_env)


# ==========================================
#           ABSEQ's DEFUALT SETTINGS
# ==========================================
DEFAULT_TOP_CLONE_VALUE = 100
RSCRIPT_PAIRING_SEPARATOR = ';'
RSCRIPT_SAMPLE_SEPARATOR = '|'
DEFAULT_MERGER = 'leehom'
WEBLOGO = 'weblogo'


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



#SAMTOOLS_PROGRAM= 'samtools'
#BGZIP = 'bgzip'
#TABIX = 'tabix'
#CLUSTALW = 'clustalw2'
