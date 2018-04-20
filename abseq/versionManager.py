import argparse
import os
import numpy
import pandas
import Bio
import datetime

from subprocess import check_output, CalledProcessError

from abseq.config import VERSION


def writeParams(args, outDir):
    """
    Writes the parameters used for analysis into analysis.params

    :param args: argparse.Namespace or dict type
            argparse namespace object, or a dict

    :param outDir: string
            output directory where analysis.params reside

    :return: string
            the filename that was produced in outDir
    """
    if isinstance(args, argparse.Namespace):
        args = vars(args)
    elif isinstance(args, dict):
        pass
    else:
        raise Exception("Unsupported parameter type {}, expecting argparse.Namespace or dict".format(type(args)))

    filename = os.path.join(outDir, "analysis.params")

    with open(filename, 'w') as out:
        out.write("AbSeq version: " + VERSION + "\n")
        out.write("IMGT version - IMGT database directory last modified time : "
                  + _get_imgt_mod_date(args['database']) + "\n")
        merger = args.get("merger", None)
        if merger:
            out.write("{} version: ".format(merger) + _get_software_version(merger) + "\n")
        out.write("IgBLAST version: " + _get_software_version('igblast') + "\n")
        out.write("pandas version: " + str(pandas.__version__) + "\n")
        out.write("numpy version: " + str(numpy.__version__) + "\n")
        out.write("biopy version: " + str(Bio.__version__) + "\n")
        out.write("FastQC version: " + _get_software_version('fastqc') + "\n")
        out.write("Clustalo version: " + _get_software_version('clustalo') + "\n")
        out.write("Executed AbSeq with the following parameters:\n")
        for key, val in args.items():
            out.write("Parameter: {:17}\tValue: {:>20}\n".format(key, str(val)))
    return os.path.basename(filename)


def _get_software_version(prog):
    """
    taken as-is from setup.py (flash version modification)
    :param prog: program name. Possible values: igblast, clustalo, fastqc, gs, leehom, flash
    :return:
    """
    try:
        if prog == 'igblast':
            retval = check_output(['igblastn', '-version']).split('\n')[1].strip().split()[2].rstrip(',')
            return str(retval)
        elif prog == 'clustalo' or prog == 'fastqc' or prog == 'gs':
            retval = check_output([prog, '--version']).strip()
            if prog == 'fastqc':
                retval = retval.split()[-1].strip().lstrip("v")
            return str(retval)
        elif prog == 'leehom':
            # leehomMulti, doesn't show version
            check_output(['which', 'leeHomMulti'])
            return "-"
        elif prog == 'flash':
            retval = check_output(['flash', '--version'])
            return (retval.split("\n")[0]).split()[-1].lstrip("v")
    except (CalledProcessError, OSError):
        return "Not found"


def _get_imgt_mod_date(fname):
    fname = os.path.abspath(os.path.expandvars(fname))
    if os.path.exists(fname):
        return str(datetime.datetime.fromtimestamp(os.path.getmtime(fname)))
    return "-"
