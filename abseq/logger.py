import logging
import os
import sys


class _Level:
    DEBUG = 'debug'
    CRIT = 'critical'
    INFO = 'info'
    WARN = 'warn'
    ERR = 'error'

    def __init__(self):
        pass


LEVEL = _Level()


def printto(stream, message, level=LEVEL.DEBUG):
    level = level.lower()

    if level not in ['debug', 'info', 'warn', 'error', 'critical']:
        raise Exception("Unknown logging level")

    if stream:
        getattr(stream, level)(message)


def setupLogger(name, logfile, flevel=logging.DEBUG, slevel=logging.WARN):
    if not os.path.exists(logfile):
        with open(logfile, 'a'): pass

    logger = logging.getLogger(name)
    logger.setLevel(flevel)

    fh = logging.FileHandler(logfile)
    fh.setLevel(flevel)

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(slevel)

    formatter = logging.Formatter("%(asctime)s (%(name)s)[%(levelname)s]: %(message)s")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

