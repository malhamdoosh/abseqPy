'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''


import os
import sys
import re
import platform


# >> helper functions
def _findWebLogo():
    """
    on unix systems, return the binary name 'weblogo'
    on windows systems, return the path to 'weblogo' script with "python" prefixed in front, EG
    "python c:\pythonN\Scripts\weblogo" (where N is 2 or 3, depending on the version of weblogo installed)

    if weblogo wasn't installed in PYTHONPATH, return 'None' (string) regardless of the OS

    :return: "weblogo" or "python c:\pythonN\Scripts\weblogo" depending on the operating system.
    Returns 'None" (string) if weblogo can't be located, regardless of the OS.
    """
    try:
        import weblogolib
        if platform.system() != 'Windows':
            return 'weblogo'
        else:
            path = weblogolib.__file__
            return "python {}".format(re.sub(r'lib.*', 'Scripts\\weblogo', path))
    except:
        return "None"


def _find_fastQC():
    """
    fastqc shebang does not work in windows, manually execute perl script using perl interpreter

    :return: "fastqc" or "perl <path>/<to>/<fastqc>" if OS is windows
    """
    if platform.system() == "Windows":
        return 'perl ' + os.path.join(os.path.expandvars("$FASTQCROOT"), 'fastqc')
    else:
        return 'fastqc'

# >> end of: helper functions


# ==========================================
#           ABSEQ's VERSION
# ==========================================

ABSEQROOT = os.path.abspath(os.path.dirname(__file__))
VERSION = '0.99.2'

# ====================================================================================
#           ABSEQ's EXTERNAL DEPENDENCIES
# ====================================================================================
# NOTE TO PROGRAMMER: IF YOU CHANGE 3rd_party TO SOME OTHER DIRECTORY NAME, MAKE SURE YOU CHANGE
# IT IN setup.py AND MANIFEST.in TOO! (just search for this comment and you'll find the exact location)
EXTERNAL_DEP_DIR = '3rd_party'
# 1. clustal omega
CLUSTALOMEGA = 'clustalo'
# 2. FASTQC
FASTQC = _find_fastQC()
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
WEBLOGO = _findWebLogo()


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
    'VH': "WGQGTXVTVSS",
    'VK': 'FGXGTKLEIK',
    'VL': 'FGXGTKLTVL'
}

FR4_CONSENSUS_DNA = {
    'VH': "TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC",
    'VK': 'TTTGGCCAGGGGACCAAGCTGGAGATCAAA',
    'VL': 'TTCGGCGGAGGGACCAAGCTGACCGTCCTA'
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
