#!env python
#Copyright (2005) Whitehead Institute for Biomedical Research
#All Rights Reserved
#Author: David Benjamin Gordon
import re, os, urllib, sys, time, getopt
import zipfile, gzip

import TAMO.paths

GLOBALS = {}


DATA = "HumanSeq Whitehead Novartis GO SGD Holstege".split()
CODE = "Clarke weblogo".split()                        


def main():
    #Main code goes here
    ARGS = getarg('args')
    datetxt   = time.strftime('%y-%m-%d',time.localtime(time.time()))
    for fcn in getarg('functions'):
        exec('%s()'%fcn)

def downloadfiles(root,urlroot,files):
    for pair in files:
        if type(pair)==type(""):
            local  = root + os.path.split(pair)[1]
            remote = urlroot + pair
            if pair.find('//')>0: remote = pair
        else:
            local  = root + pair[1]
            remote = urlroot + pair[0]
            if pair[0].find('//')>0: remote = pair[0]
        download2local(remote,local)
    
def HumanSeqData():
    root    = TAMO.paths.HumanSeqdir
    urlroot = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg17/bigZips/'

    remote = urlroot + 'chromFa.zip'
    local =  root    + 'chromFa.zip'

    download2local(remote,local)

    print "\nUNZIPPING %s"%local
    sys.stdout.flush()
    try:
        code = os.system('cd %s && unzip %s'%(root,local))
    except:
        code = 1
    if code:
        print "\nError unzipping:  attempted to execute command:"
        print 'cd %s && unzip %s'%(root,local)
        sys.exit(1)

def HolstegeData():
    root    = TAMO.paths.Holstegedir
    urlroot = 'http://web.wi.mit.edu/young/pub/'
    files   = ['orf_transcriptome.txt']
    downloadfiles(root,urlroot,files)

def WhiteheadData():
    root    = TAMO.paths.Whiteheaddir
    urlroot = 'http://fraenkel.mit.edu/Harbison/release_v24/'
    files   = ['Harbison_Gordon_yeast_v9.11.csv.gz',
               'Yeast6kArray.tgz',
               'Human13kArray.tgz']

    downloadfiles(root,urlroot,files)
    for file in files:
        gunzip(root,file)

def NovartisData():
    root    = TAMO.paths.Novartisdir
    urlroot = 'http://wombat.gnf.org/downloads/'
    files = ['GNF1Hdata.zip',
             'gnf1h-anntable.txt.gz']
    downloadfiles(root,urlroot,files)

    unzip(root,files[0],lambda x: x.replace('\t',','))


def GOData():
    root    = TAMO.paths.SGDdir
    urlroot = 'ftp://genome-ftp.stanford.edu/pub/yeast/data_download/' 
    files = ['literature_curation/go_slim_mapping.tab']
    downloadfiles(root,urlroot,files)

def weblogoCode():
    root = TAMO.paths.weblogodir 
    urlroot = 'http://weblogo.berkeley.edu/'
    files   = ['release/weblogo.2.8.1.tar.gz', 'LICENSE']
    downloadfiles(root,urlroot,files)
    gunzip(root,'weblogo.2.8.1.tar.gz')
    print "WEBLOGO -- DOWNLOADED, BUT REQUIRES MANUAL INSTALL"
    print "GO TO %S DIRECTORY FOR INSTRUCTIONS"%root
    

def ClarkeCode():
    root    = TAMO.paths.TAMOroot + 'util/Clarke/'
    urlroot = 'ftp://ftp.bs.jhmi.edu/users/nclarke/MNCP/'
    files   = ['ROC_AUC.py',
               'multi_test.py',
               'narke_test.py',
               'readme.txt']
    downloadfiles(root,urlroot,files)
            

def SGDData():
    root    = TAMO.paths.SGDdir
    urlroot = 'ftp://genome-ftp.stanford.edu/pub/yeast/data_download/' 
    files = ['chromosomal_feature/SGD_features.tab',
             'chromosomal_feature/dbxref.tab',
             'chromosomal_feature/chromosome_length.tab',
             'sequence/GenBank/yeast_nrpep.fasta.gz',
             'sequence/genomic_sequence/orf_protein/orf_trans_all.fasta.gz',
             ('http://yeastgfp.ucsf.edu/allOrfData.txt','Huh_Nature_2003.tab')
             ]

    chrs = '01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 mt'.split()

    files.extend( ['sequence/NCBI_genome_source/chr%s.fsa'%x for x in chrs] )

    downloadfiles(root,urlroot,files)

    from TAMO.seq import Fasta
    
    print "Assembling yeast genome sequence files into a single file (NCBI_yeast_genome.fsa)"
    D = {}
    for chr in chrs:
        _d = Fasta.load('%s/chr%s.fsa'%(TAMO.paths.SGDdir,chr))
        id, seq = _d.items()[0]
        if chr[0] == '0': chr = chr[1]
        D['chr%s  %s'%(chr,id)] = seq
    Fasta.write(D, TAMO.paths.SGDdir + 'NCBI_yeast_genome.fsa')

def gunzip(root,file):
    if file.find('.tgz') > 0:
        try:
            code = os.system('cd %s && tar xzvf %s'%(root,file))
        except:
            code = 1
        if code:
            print "\nError un-tar-gzipping:  Attempted command:"
            print 'cd %s && tar xzvf %s'%(root,file)
    else:
        newfile = file.replace('.gz','')
        print " gnu-unzipping %s/%s -> %s/%s"%(root,file,root,newfile)
        F = gzip.open(root+file)
        O = open(root+newfile,'w')
        O.write(F.read())
        F.close()
        O.close()
    
    

def unzip(root,file,FCN=None):
    print " Unzipping %s/%s .... "%(root,file)
    F = zipfile.ZipFile(root+file)
    filenames = F.namelist()
    for filename in filenames:
        print '\t',filename
        txt = F.read(filename)
        if FCN: txt = FCN(txt)
        path,tail = os.path.split(filename)
        if path and not os.path.exists(root+path): os.makedirs(root+path)
        O = open(root+filename,'w')
        O.write(txt)
        O.close()

def download2local(url,localfile):
    if os.path.exists(localfile):
        print "!!! The following file already exists.  Delete manually if ",
        print "    you wish to overwrite:\n!!! %s\n"%localfile
        return
    path,tail = os.path.split(localfile)
    if not os.path.exists(path):  os.makedirs(path)
    print "\n  Downloading from:   %s"%url
    print   "                to:   %s ..."%localfile,
    sys.stdout.flush()
    try:
        urllib.urlretrieve(url,localfile, dashes)
        print
    except:
        print 
        errstr = "#! Could not download \n\nfrom: %s\n  to: %s.  \n\n"%(url,localfile)
        errstr +="#! Check URL for correctness, and path for permissions\n"
        print errstr
        urllib.urlretrieve(url,localfile, dashes)
        usage()

def dashes(numblocks,blocksize,totsize):
    total_k   = totsize/1024.0
    current_k = numblocks*blocksize/1024.0
    if current_k > total_k: current_k = total_k

    if numblocks == 0: print "(Total: %8dk) : "%(total_k),
    else:              sys.stdout.write('\010'*9)
    sys.stdout.write('%8dk'%(current_k))
    
    #curp20  = float(numblocks)   * float(blocksize) / float(totsize) * 20
    #lastp20 = float(numblocks-1) * float(blocksize) / float(totsize) * 20
    #if int(lastp20) != int(curp20):
    #    sys.stdout.write('-')
    #    sys.stdout.flush()

def download2lines(url):
    try:
        F = orllib.URLopener().open(url)
        lines = F.readlines()
        F.close()
        return F
    except:
        usage("#! Could not download %s.  \n#! Check URL for correctness\n");

def usage(txt=''):
    if txt: print "Error: %s"%txt
    print 'Usage: %s [--all  | --<data/code source>]'%(
        re.sub('^.*/','',sys.argv[0]))
    print ''
    print '   Data sources include:'
    for d in DATA:    print '         --%s'%d
    print '   Code sources include:'
    for c in CODE:    print '         --%s'%c
    print '   Note, some datasets (e.g. HumanSeq) are very large, and may take a long time to download'
    #print '           [-w window (100)]'
    sys.exit(1)

def parse_opts():
    global GLOBALS
    short_opts = ''
    long_opts  = ['all']
    long_opts.extend(DATA[:])
    long_opts.extend(CODE[:])
    
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()

    dDATA = ['--%s'%x for x in DATA]
    dCODE = ['--%s'%x for x in CODE]

    GLOBALS['args'] = args
    GLOBALS['file'] = ''
    GLOBALS['functions'] = []
    for opt,value in opts:
        if opt in dDATA:                  GLOBALS['functions'].append('%sData'%opt.replace('--',''))
        if opt in dCODE:                  GLOBALS['functions'].append('%sCode'%opt.replace('--',''))
        if opt == '--all':                GLOBALS['all']   = 1

    if getarg('all'):
        GLOBALS['functions']    =   ['%sData'%x for x in DATA]
        GLOBALS['functions'].extend(['%sCode'%x for x in CODE])
        

def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None

if __name__ == '__main__': 
    parse_opts()
    print "#" + ' '.join([x.replace(' ','\ ') for x in sys.argv])
    main()                                  


