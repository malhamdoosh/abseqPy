'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''


import os
import sys
import platform

# ==========================================
#           ABSEQ's VERSION
# ==========================================
ABSEQROOT = os.path.abspath(os.path.dirname(__file__))
VERSION = '1.1.15'

# ====================================================================================
#           ABSEQ's EXTERNAL DEPENDENCIES
# ====================================================================================
# NOTE TO PROGRAMMER: IF YOU CHANGE 3rd_party TO SOME OTHER DIRECTORY NAME, MAKE SURE YOU CHANGE
# IT IN setup.py AND MANIFEST.in TOO! (just search for this comment and you'll find the exact location)
EXTERNAL_DEP_DIR = '3rd_party'
# 1. clustal omega
CLUSTALOMEGA = 'clustalo'
# 2. FASTQC
FASTQC = 'fastqc'
# 3. mergers
LEEHOM = 'leeHomMulti'
FLASH = 'flash'
PEAR = 'pear'
# 4.a IgBlastn
IGBLASTN = 'igblastn'
# 4.b IgBlastp
IGBLASTP = 'igblastp'




# ==========================================
#           ABSEQ's DEFUALT SETTINGS
# ==========================================
DEFAULT_TOP_CLONE_VALUE = 'inf'
DEFAULT_MERGER = 'leehom' if platform.system() != "Windows" else 'flash'
DEFAULT_TASK = 'abundance'
WEBLOGO = 'weblogo'


#  === Chain type: K ===
# VDVRPRDQGGNQ                   WTFGQGTKVEIK                   GRSAKGPRWKSN
# CTLLARGPSWRSN                  VHFWPGDQAGDQ                   YTFGQGTKLEIK
# VHFWPGDQAGDQ                   CTFGQGTKLEIK @                 ALLARGPSWRSN
# CTVLARGPSWRSN                  VQFWPGDQAGDQ                   YSFGQGTKLEIK
# CAVLARGPSWRSN                  VQFWPGDQAGDQ                   CSFGQGTKLEIK
# IHFRPWDQSGYQ                   FTFGPGTKVDIK                   SLSALGPKWISN
# AHFRRRDQGGDQ                   LTFGGGTKVEIK                   SLSAEGPRWRSN
# AHVRRRDQGGDQ                   LTFGGGTKVEIK                   SRSAEGPRWRSN
# DHLRPRDTTGD*                   ITFGQGTRLEIK                   SPSAKGHDWRLN
#  === Chain type: L ===
# LCLRNWDQGHRP                   YVFGTGTKVTVL                   MSSELGPRSPS*
# CGIRRRDQADRP                   VVFGGGTKLTVL @                 WYSAEGPS*PS*
# CGIRRRDQADRP                   VVFGGGTKLTVL                   WYSAEGPS*PS*
# LGVRRRDQADRP                   WVFGGGTKLTVL                   GCSAEGPS*PS*
# FCIWWRNPADHF                   FVFGGGTQLIIL                   LYLVEEPS*SF*
# LGVW*GDRADRP                   WVFGEGTELTVL                   GCLVRGPS*PS*
# LGVW*GDGADRP                   WVFGEGTELTVL                   GCLVRGRS*PS*
# *CVRQWHQGDRP                   NVFGSGTKVTVL                   MCSAVAPR*PSS
# CCVRRRHPADRP                   AVFGGGTQLTVL                   LCSEEAPS*PSS
# CCVRRRHPADRP                   AVFGGGTQLTAL                   LCSEEAPS*PPS
# consensus protein of HV http://discovery.ucl.ac.uk/15808/1/15808.pdf
FR4_CONSENSUS = {
    'hv': "WGQGTXVTVSS",
    'kv': 'FGXGTKLEIK',
    'lv': 'FGXGTKLTVL'
}

FR4_CONSENSUS_DNA = {
    'hv': "TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC",
    'kv': 'TTTGGCCAGGGGACCAAGCTGGAGATCAAA',
    'lv': 'TTCGGCGGAGGGACCAAGCTGACCGTCCTA'
}

# directory naming
AUX_FOLDER = 'auxiliary'
HDF_FOLDER = 'hdf'


GB = (1024.**3)
# sorry darwin people, you need psutil because sysconf can't locate 'sc_phys_pages'
if sys.platform == 'darwin' or platform.system() == "Windows":
    from psutil import virtual_memory
    mem = virtual_memory()
    MEM_GB = mem.total/GB
else:
    tmp = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
    MEM_GB = tmp/GB


# SAMTOOLS_PROGRAM= 'samtools'
# BGZIP = 'bgzip'
# TABIX = 'tabix'
# CLUSTALW = 'clustalw2'
