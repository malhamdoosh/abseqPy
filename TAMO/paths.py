'''
TAMO.paths: Path information for TAMO Bioninformatics Modules

   NOTE: When adding entries, append a "/" to directory names
   to facilitate TAMO.paths.variable + "filename" convention

   Instructions for downloading Motif Discovery Programs
   (AlignACE, MEME, & MDscan) can be found near the end
   of the file.

Copyright (2005) Whitehead Institute for Biomedical Research
All Rights Reserved
Author: David Benjamin Gordon

'''
import os,sys
import TAMO.localpaths
TAMOroot = TAMO.localpaths.TAMOroot
TAMOdata = TAMO.localpaths.TAMOdata

#TAMOroot  = '/etc/TAMO/'      #Installation directory for TAMO (might be /usr/lib/python2.4/site-packages/TAMO/)
#TAMOdata  = '/etc/TAMOdata/'  #Directory for storing data downloaded over the internet

THEMEroot = TAMOroot

########################################################################
#                      Paths for TAMO modules                          #
########################################################################

FSAdir  = TAMOdata + 'fsafiles/'
BGdir = TAMOdata + 'fsafiles/'

########################################################################
#              Paths for data from the Whitehead Institute             #
#                  Largely from Fraenkel and Young Labs                #
########################################################################

Whiteheaddir    = TAMOdata + 'Whitehead/'
Yeast6kArraydir  = Whiteheaddir + 'Yeast6kArray/'
Human13kArraydir = Whiteheaddir + 'Human13kArray/'

########################################################################
#                      Paths for SGD module                            #
#                      (also used by GO.py)                            #
########################################################################

SGDdir = TAMOdata + 'SGD/'

########################################################################
#                Paths for Human Sequence Information                  #
########################################################################

HumanSeqdir = TAMOdata + 'HumanSeq/'

########################################################################
#                    Paths for Novartis Module                         #
########################################################################

Novartisdir = TAMOdata + 'Novartis/'


########################################################################
#                    Paths for Holstege Data                           #
########################################################################

Holstegedir = TAMOdata + 'Holstege/'


########################################################################
#     Paths (and instructions for other Motif Discovery Programs       #
########################################################################


'''
# AlignACE #

AlignACE can be downloaded after a click-through licence agreement at:

http://atlas.med.harvard.edu/download/

NOTE: It _does_ work to download the "linux" AlignACE package and compile under cygwin

Replace the following line with the name of the directory into which you extracted (and built?) AlignACE below:
'''

AlignACEdir = TAMOdata + 'MDprogs/alignace2004/'

'''
# MEME #

MEME can be downloaded from

http://meme.sdsc.edu/meme/website/meme-download.html

NOTE: It _does_ work to build MEME under cygwin.

Replace the following line with the name of the directory into which you extracted (and built?) MEME below:
'''

MEMEdir = TAMOdata + 'MDprogs/meme.3.0.13/'

'''
# MDscan #

MDscan is available only for linux (as far as I know) after faxing a licence agreement
to the Brutlag Group.  Instructions are at:

http://motif.stanford.edu/distributions/

Replace the following line with the name of the directory into which you extracted MDscan below:
'''

MDscandir = TAMOdata + 'MDprogs/MDscan/'

########################################################################
#     Weblogo                                                          #
########################################################################

"""
If you've installed weblogo put the path to the executable here.  Weblogo/Seqlogo
can be downloaded from http://weblogo.berkeley.edu/.
"""

weblogodir = TAMOdata + 'weblogo/'




def CHECK(filelist,arg,note=''):
    """
    Check whether the files can be found in the expected locations.  If not,
    suggest how they might be retrieved.
    """
    if type(filelist) != type([]): filelist = [filelist]
    for file in filelist:
        if note:
            exceptiontxt =  '\n   CANNOT FIND FILE %s\n   EXAMINE TAMO/paths.py '%(file)
            exceptiontxt += 'FOR INFORMATION REGARDING %s'%(note)
            exceptiontxt += '   (If file does exist, check permissions on file and directory)'
        else:
            exceptiontxt  = '   CANNOT FIND FILE %s\nINVOKE:\n\t    GetDataFiles.py --%s'%(file,arg)
            exceptiontxt += '\n   (If file does exist, check permissions on file and directory)'

        if not os.path.exists(file):
            raise Exception, exceptiontxt
        

