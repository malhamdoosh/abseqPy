'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import print_function
import os
import sys
import argparse
from config import VERSION
from numpy import Inf
from os.path import abspath

def extractRanges(strRanges, expNoRanges=2):
    """
    Returns the range given an input
    :param strRanges: string range, allowed format: 0-10,23-25 or 0-10 or 0
    :param expNoRanges: maximum number of allowed ranges, eg: expNoRange=1 implies only 1 range, expNoRange=2
    implies 0-10,20-25 is allowed, and so on
    :return: a nested list of ranges
    """
    numRanges = []
    ranges = strRanges.split(',')
    if (len(ranges) > expNoRanges):
        raise Exception("Number of bitScore, alignLen and sstart ranges should match the number of files")
        
    for i in range(len(ranges)):
        scores = ranges[i].split('-')
        if (len(scores) == 2):
            numRanges.append([float(scores[0]), float(scores[1])])
            if (numRanges[-1][0] >= numRanges[-1][1]):
                raise Exception("Invalid ranges " + strRanges)
                
        else:
            numRanges.append([float(scores[0]), Inf])

    # BUGSQ: if (len(numRanges) == 1 < expNoRanges):
    if len(numRanges) < expNoRanges:
        numRanges = numRanges * expNoRanges
    return numRanges

PROGRAM_VALID_ARGS = ['-task', '-chain', '-name',
                 '-f1', '-f2', '-fmt', '-o', '-merge', '-merger',  
                 '-seqtype', '-threads', '-db',  
                  '-bitscore', '-alignlen', '-sstart', '-actualqstart',     
                  '-trim5' ,'-trim3',  '-fr4cut',            
                   '-sites',
                  '-primer', 
                  '-5end', '-3end', 
                  '-5endoffset',
                 '-upstream',
                 '-report_interim'
                  ]

# def parseArgs(args):
#     argVals = {}
#     argVals["cmd"] = ' '.join(args)
#     for i in range(1, len(args), 2):
# #         print(args[i].lower(), PROGRAM_VALID_ARGS)
#         if (args[i].lower() in PROGRAM_VALID_ARGS):
# #             print(args[i], args[i+1])
#             argVals[args[i][1:].lower()] = args[i+1]
#         else:
#             print(args[i], PROGRAM_VALID_ARGS)
#             raise Exception(args[i] + ' is invalid argument.')
#
#     if (argVals.get('chain', None) is None):
#         argVals['chain'] = "hv"
#     else:
#         argVals['chain'] = argVals['chain'].lower()
#     if (argVals.get('task', None) is None):
#         argVals['task'] = 'abundance'
#     else:
#         argVals['task'] = argVals['task'].lower()
#
#     if (argVals.get('seqtype', None) is None):
#         argVals['seqtype'] = 'dna'
#
#     # input sequencing files
#     if (argVals.get('f1', None) is None):
#         raise Exception("One sequence file at least must be provided.")
#     elif not os.path.exists(argVals['f1']):
#         raise Exception("-f1 file not found!")
#     else:
#         argVals['f1'] = abspath(argVals['f1'])
#
#     if (argVals.get('f2', None) is None):
#         argVals['f2'] = None
#     elif not os.path.exists(argVals['f2']):
#         raise Exception("-f2 file not found!")
#     else:
#         argVals['f2'] = abspath(argVals['f2'])
#
#     # merging options
#     if (argVals.get('merge', None) is None):
#         argVals['merge'] = 'no'
#     else:
#         argVals['merge'] = argVals['merge'].lower()
#     if (argVals['merge'] == 'yes' and (argVals['f2'] is None)):
#         raise Exception("The merger requires two sequence files (use both -f1 and -f2).")
#     else:
#         if (argVals.get('merger', None) is None):
#             argVals['merger'] = 'flash'
#         else:
#             argVals['merger'] = argVals['merger'].lower()
#
#     ## output directory and parse sample name
#     if (argVals.get('o', None) is None):
#         argVals['o'] = './'
#     if (argVals.get('name', None) is not None):
#         argVals['o'] += "/" + argVals['name']
#     else:
#         f1name = argVals['f1'].split('/')[-1]
#         if (f1name.find("_R") != -1 and
#             (argVals['merge'] == 'yes' or argVals['task'] == "fastqc")):
#             ext = '_' + f1name.split('_')[-1]
#         else:
#             ext = f1name[f1name.find('.'):]
#         argVals['o'] += "/"+f1name.replace(ext, '')
#         sampleName = f1name.split("_")[0] + '_'
#         sampleName += f1name.split("_")[-1].split('.')[0]
#         argVals['name'] = sampleName
#     argVals['o'] = (abspath(argVals['o']) + '/').replace('//', '/')
#     os.system("mkdir -p " + argVals['o'])
#
#
#
#     # (argVals['merge'] == 'no' and noFiles == 2) or
#
#
#     if (argVals.get('fmt', None) is None):
#         argVals['fmt'] = 'fastq'
#     else:
#         argVals['fmt'] = argVals['fmt'].lower()
#
#     if (argVals.get('bitscore', None) is None):
#         argVals['bitscore'] = [0, Inf]
#     else:
#         argVals['bitscore'] = extractRanges(argVals['bitscore'])[0]
#
#     if (argVals['task'] in ['secretion', '5utr']):
#         if (argVals.get('upstream', None) is None):
#             argVals['upstream'] = [1, Inf]
#         else:
#             argVals['upstream'] = extractRanges(argVals['upstream'], 1)[0]
#     if (argVals['task'] in ['rsa', 'rsasimple']):
#         if (argVals.get('sites', None) is None):
#             raise Exception("Restriction sites should be provided.")
#         argVals['sites'] = abspath(argVals['sites'])
#     if (argVals['task'] in ['diversity', 'productivity', 'all']):
#         # actualqstart is a 1-based index
#         if argVals.get('actualqstart', None) is not None:
#             argVals['actualqstart'] = int(argVals['actualqstart']) - 1
#         else:
#             argVals['actualqstart'] = -1
#         if (argVals.get('fr4cut', None) is None):
#             argVals['fr4cut'] = False
#         else:
#             argVals['fr4cut'] = (argVals['fr4cut'].upper() == 'YES')
#     if (argVals.get('trim5', None) is not None):
#         argVals['trim5'] = int(argVals['trim5']) - 1
#     else:
#         argVals['trim5'] = 0
#     if (argVals.get('trim3', None) is not None):
#         argVals['trim3'] = int(argVals['trim3'])
#     else:
#         argVals['trim3'] = 0
#     if (argVals['task'] == 'primer'):
#         argVals['5end'] = argVals.get('5end', None)
#         if argVals['5end']:
#             argVals['5end'] = abspath(argVals['5end'])
#         argVals['3end'] = argVals.get('3end', None)
#         if argVals['3end']:
#             argVals['3end'] = abspath(argVals['3end'])
#         if (argVals.get('5endoffset', None) is None):
#             argVals['5endoffset'] = 0
#         else:
#             argVals['5endoffset'] = int(argVals['5endoffset'])
#
#     if (argVals.get('sstart', None) is None):
#         argVals['sstart'] = [1, Inf]
#     else:
#         argVals['sstart'] = extractRanges(argVals['sstart'])[0]
#
#     if (argVals.get('alignlen', None) is None):
#         argVals['alignlen'] = [0, Inf]
#     else:
#         argVals['alignlen'] = extractRanges(argVals['alignlen'])[0]
#
#     if (argVals.get('primer', None) is None):
#         argVals['primer'] = -1
#     else:
#         argVals['primer'] = int(argVals['primer'])
#
#     if (argVals.get('db', None) is None):
#         argVals['db'] = '$IGBLASTDB'
#     else:
#         argVals['db'] = os.path.abspath(argVals['db'])
#
#     if (argVals.get('threads', None) is None):
#         argVals['threads'] = 8
#     else:
#         argVals['threads'] = int(argVals['threads'])
#
#     if (argVals.get('report_interim', None) is None):
#         argVals['report_interim'] = False
#     else:
#         argVals['report_interim'] = (argVals['report_interim'].upper() == 'YES')
#     argVals['log'] = argVals['o'] + argVals['name'] + ".log"
#
#     return argVals
#


def printUsage(parser, additional_msg=None):
    parser.print_help()
    if additional_msg is not None:
        print(additional_msg, file=sys.stderr)
    sys.exit(0)


def parseArgs():
    """
    Parses sys.argv's arguments and sanitize them. Checks the logic of arguments so calling program does not have
    to do any logic checking after this call.
    :return: argparse namespace object, using dot notation to retrieve value: args.value
    """
    parser, args = parseCommandLineArguments()

    # canonicalize all values
    args.task = args.task.lower()
    args.seqtype = args.seqtype.lower()
    args.chain = args.chain.lower()
    args.fmt = args.fmt.lower()
    args.merger = args.merger.lower()

    # check for f1, f2 file existence and expand path
    if not os.path.exists(args.f1):
        raise Exception("-f1 file not found!")
    else:
        args.f1 = abspath(args.f1)
    if args.f2 is not None and not os.path.exists(args.f2):
        raise Exception("-f2 file not found!")
    elif args.f2 is not None:
        args.f2 = abspath(args.f2)

    # check logic between f1, f2 and merger, setting default merger to flash
    if args.merger is not None and args.f2 is None:
        raise Exception("The merger requires two sequence files (use both -f1 and -f2 option)")
    if args.merger is None and args.f2 is not None:
        # flash is default merger, as per help message
        args.merger = 'flash'

    # appending analysis name to output directory filenames
    if args.name is not None:
        args.outdir += ("/" + args.name)
    else:
        f1name = args.f1.split("/")[-1]
        if f1name.find("_R") != -1 and (args.merger is not None or args.task.lower() == "fastqc"):
            ext = '_' + f1name.split("_")[-1]
        else:
            ext = f1name[f1name.find("."):]
        args.outdir += "/" + f1name.replace(ext, "")
        sampleName = f1name.split("_")[0] + "_"
        sampleName += f1name.split("_")[-1].split(".")[0]
        args.name = sampleName
    args.outdir = (abspath(args.outdir) + "/").replace("//", "/")

    # silently ignore creation of out directory if already exists
    if not os.path.exists(args.outdir):
        os.mkdirs(args.outdir)

    # setting default value of bitscores if not provided, else extract the string ranges provided
    args.bitscore = [0, Inf] if args.bitscore is None else extractRanges(args.bitscore)[0]

    # setting default values for upstream
    if args.task in ['secretion', '5utr']:
        args.upstream = [1, Inf] if args.upstream is None else extractRanges(args.upstream, 1)[0]

    # confirm that file to sites is provided
    if args.task in ['rsa', 'rsasimple']:
        if args.sites is None:
            printUsage(parser, "Restriction sites should be provided if --task rsa or --task rsasimple")
            # exit from print_usage, never reaching here
        args.sites = abspath(args.sites)

    # provided actualqstart is converted to 0-base from 1-based index, -1 is checked later on for default value
    if args.task in ['diversity', 'productivity', 'all']:
        args.actualqstart = (args.actualqstart - 1) if args.actualqstart is not None else -1

    # BUGSQ: if user provided value = 0, what happens?: here, only subtract 1 if args.trim5 isn't default 0, or if user
    # didn't provide 0, since the other file that uses this parameter didn't check for negative values
    args.trim5 -= args.trim5 > 0            # again, transform args.trim5 to 0-based if provided value is > 0, else 0

    # retrieve filenames for primer analysis on 5' and 3' end
    if args.task == 'primer':
        args.primer5end = abspath(args.primer5end) if args.primer5end is not None else None
        args.primer3end = abspath(args.primer3end) if args.primer3end is not None else None

    args.sstart = [1, Inf] if args.sstart is None else extractRanges(args.sstart)[0]
    args.alignlen = [0, Inf] if args.alignlen is None else extractRanges(args.alignlen)[0]
    args.database = abspath(args.database) if args.database is not None else "$IGBLASTDB"
    args.log = args.outdir + args.name + ".log"

    return args


def parseCommandLineArguments():
    """
    parses commandline arguments for AbSeq
    :param argv: sys.argv
    :return: parser object, can be indexed for flag values
    """
    parser = argparse.ArgumentParser(description='AbSeq antibody library sequencing quality control pipeline',
                                     prog="AbSeq")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument('-c', '--chain', default="hv", help="Chain type [default=hv]",
                        choices=['hv', 'lv', 'kv'])
    parser.add_argument('-t', '--task', default="abundance", help="Analysis task, supported tasks: \
                                                    all, annotate, abundance, \
                                                    diversity, fastqc, productivity, primer, 5utr, rsasimple, rsa, \
                                                    seqlen, secretion, seqlenclass [default=abundance]",
                                                    choices=["all", "annotate", "abundance", "diversity",
                                                             "fastqc", "productivity", "primer", "5utr", "rsasimple",
                                                             "rsa", "seqlen", "secretion", "seqlenclass"])
    parser.add_argument('-s', '--seqtype', default='dna', help="Sequence type, supported seq type: dna or protein \
                                                                    [default=dna]",
                                                        choices=["dna", "protein"])
    parser.add_argument('-f1', '--file1', dest="f1", required=True, help="Fully qualified path to sequence file 1")
    parser.add_argument('-f2', '--file2', dest="f2", help="Fully qualified path to sequence file 2,"
                                                          " leave out if sequences are not paired end", default=None)
    # in parseArgs, we change None to flash by default if there's a -f2 option
    parser.add_argument('-m', '--merger', help="Choice between different mergers. Omit this if no -f2 option"
                                                " is specified [default=flash]",
                        default=None,
                        choices=['leehom', 'flash', 'pear'])
    parser.add_argument('-o', '--outdir', help="Output directory [default = current working directory]", default="./")
    parser.add_argument('-n', '--name', help="Name of analysis [default = name of Sequence file 1/2]", default=None)
    parser.add_argument('-f', '--format', dest="fmt", help="Format of input file",
                        default="fastq", choices=['fasta', 'fastq'])
    parser.add_argument('-b', '--bitscore', help="Bitscore threshold to apply on V gene, accepted format: num1-num2"
                                                 " [default=[0, inf)]", default=None)
    parser.add_argument('-t5', '--trim5', help="Number of nucleotides to trim on the 5'end of V gene [default=0]",
                        default=0, type=int)
    parser.add_argument('-t3', '--trim3', help="Number of nucleotides to trim on the 3'end of V gene [default=0]",
                        default=0, type=int)
    parser.add_argument('-ss', '--sstart', help="Filter sequences that do not fall into this start range, accepted"
                                                " format: num1-num2 [default=[1, inf)]", default=None)
    parser.add_argument('-al', '--alignlen', help="Filter sequences that do not fall into this alignment length range,"
                                                  " accepted format: num1-num2 [default=[0, inf)]", default=None)
    parser.add_argument('-p', '--primer', help="Not implemented yet [default=-1]", default=-1, type=int)
    parser.add_argument('-d', '--database', help="Specify fully qualified path to germline database "
                                                           "[default=$IGBLASTDB], type echo $IGBLASTDB in command line"
                                                            " to see your default database used by AbSeq",
                        default=None)
    parser.add_argument('-q', '--threads', help="Number of threads to use (spawns separate processes)[default=8]",
                        type=int, default=8)
    parser.add_argument('-r', '--report-interim', help="Not implemented yet [default= no report]",
                        default=False, type=bool)
    parser.add_argument('-u', '--upstream', help="Range of upstream sequences, secretion signal analysis and 5UTR"
                                                 " analysis [default=[1, inf)]", default=None)
    parser.add_argument('-st', '--sites', help="Fully qualified pathname to restriction sites file, required if"
                                               " --task rsa or --task rsasimple is specified", default=None)
    parser.add_argument('-qs', '--qstart', dest="actualqstart", help="Specify starting position of query V gene during"
                                                                     " alignment (1-based indexing) [default=1]",
                        default=None, type=int)
    parser.add_argument('-f4c', '--fr4cut', help="Cut framework 4 region [default = no cuts]",
                        type=bool, default=False)
    parser.add_argument('-p3', '--primer3end', help="Fully qualified path to primer 3' end file", default=None)
    parser.add_argument('-p5', '--primer5end', help="Fully qualified path to primer 5' end file", default=None)
    parser.add_argument('-p5off', '--primer5endoffset', help="Number of nucleotides for 5' end offset [default=0]",
                        default=0, type=int)
    return parser, parser.parse_args()
