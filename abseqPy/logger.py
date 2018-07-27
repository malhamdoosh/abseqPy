import logging
import os
import sys


class _Level:
    DEBUG = 'debug'
    CRIT = 'critical'
    INFO = 'info'
    WARN = 'warn'
    ERR = 'error'
    EXCEPT = 'exception'

    def __init__(self, streamLevel=logging.WARN, fileLevel=logging.DEBUG):
        self.streamLevel = streamLevel
        self.fileLevel = fileLevel


# LOG FILE HEADER - displayed on the top of each *.log file depending on analysis task
_BANNER = {
    'all': "Running the complete QC pipeline",
    'fastqc': "Sequencing QC Analysis",
    'annotate': "Clone Identification and Classification",
    'abundance': "IGV Abundance and QC Plots",
    'productivity': "Clone Productivity Analysis",
    'diversity': "Diversity Analysis",
    'secretion': "Secretion signal Analysis",
    '5utr': "5'UTR analysis",
    'rsasimple': "Simple Restriction Sites Analysis",
    'rsa': "Comprehensive Restriction Sites Analysis",
    'primer': "Primer Specificity Analysis",
    'seqlen': "Sequence Length Distribution",
    'default': "Running AbSeq on all sample(s)"
}

LEVEL = _Level()


def printto(stream, message, level=LEVEL.DEBUG):
    level = level.lower()

    if level not in [_Level.DEBUG, _Level.CRIT, _Level.INFO, _Level.WARN, _Level.ERR, _Level.EXCEPT]:
        raise Exception("Unknown logging level")

    if stream:
        getattr(stream, level)(message)


def setupLogger(name, task, logfile, stream=sys.stdout, flevel=LEVEL.fileLevel, slevel=LEVEL.streamLevel):
    with open(logfile, 'a') as fp:
        fp.write(formattedTitle(task) + '\n')

    datetimefmt = "%Y-%m-%d %H:%M:%S"

    logger = logging.getLogger(name)
    logger.setLevel(flevel)

    fh = logging.FileHandler(logfile)
    fh.setLevel(flevel)

    ch = logging.StreamHandler(stream=stream)
    ch.setLevel(slevel)

    formatter = logging.Formatter("%(asctime)s (%(name)s)[%(levelname).4s]: %(message)s", datefmt=datetimefmt)
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)


def formattedTitle(task, defaultTitle=False):
    if defaultTitle:
        title = _BANNER['default']
    else:
        title = _BANNER.get(task, None)
    if title is None:
        raise Exception("Unknown task requested. Available tasks are: {}".format(','.join(_BANNER.keys())))
    string = "-" * 100 + '\n'
    string += "|" + " " * 98 + "|\n"
    string += "|" + " " * ((98 - len(title)) / 2) + title + " " * (
            (98 - len(title)) / 2 + (98 - len(title)) % 2) + "|\n"
    string += "|" + " " * 98 + "|\n"
    string += "-" * 100 + '\n'
    return string
