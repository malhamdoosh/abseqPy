import os
import sys
import subprocess
import shlex

from abseqPy.config import ABSEQROOT, EXTERNAL_DEP_DIR
from abseqPy.config import MEM_GB


# temporarily overrides PATH variable with EXTERNAL_DEP_DIR/bin, IGBLASTDB and IGDATA (if they exist)
class PriorityPath:
    def __init__(self):
        self.updated = False
        self.old_env = os.environ.copy()
        _env = os.environ.copy()

        # if the BIN directory exists, append it to the front of PATH variable
        override_path = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'bin')) + os.path.sep
        if os.path.exists(override_path):
            _env['PATH'] = override_path + os.pathsep + _env['PATH']
            self.updated = True

        # if the igdata dir exists, override it irrespective of if there's already a IGDATA env
        override_igdata = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'igdata')) + os.path.sep
        if os.path.exists(override_igdata):
            _env['IGDATA'] = override_igdata
            self.updated = True

        # if the igdb dir exists, override it irrespective of if there's already a IGBLASTDB env
        override_igdb = os.path.abspath(os.path.join(ABSEQROOT, EXTERNAL_DEP_DIR, 'databases')) + os.path.sep
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


def hasLargeMem(size=16):
    """
    tells if system has memory strictly larger than specified size in GB

    :param size: unit GB
    :return: bool. true if virtual_memory > size
    """
    return MEM_GB > size


class CommandLine:
    def __init__(self, exe, *args, **kwargs):
        self._exe = exe
        self._kwargs = kwargs
        self._args = args
        self._ext = ""

    def append(self, string):
        """
        appends extra arguments behind the string-ified command. This is
        useful for edge cases that aren't covered by this class.
        For example, to mix long and short options, to add "=" in options.

        :param string: string. Will be appended directly behind the command when executed
        :return: self

        >>> cmd = ShortOpts("ls", "$HOME", l="").append("-a").append("-R")
        >>> str(cmd)
        'ls $HOME -l -a -R'
        >>> cmd.append("-f")
        ls $HOME -l -a -R -f
        """
        self._ext = (self._ext + " " + string)
        # allow chaining
        return self

    def __call__(self, stdout=sys.stdout, stderr=sys.stderr):
        """
        executes the built command. Raises subprocess.CalledProcessError if
        something goes wrong

        :param stdout: std output stream. A value of None will flush it to /dev/null
        :param stderr: std error stream. A value of None will flush it to /dev/null
        :return: None

        >>> # execute python -m this
        >>> cmd = ShortOpts("python", m='this')
        >>> cmd
        python -m this
        >>> tmpdir = getfixture("tmpdir")
        >>> with tmpdir.join("tmp.txt").open("w") as fp:
        ...    cmd(stdout=fp)
        >>> with tmpdir.join("tmp.txt").open("r") as fp:
        ...    fp.readlines()[0].strip()
        'The Zen of Python, by Tim Peters'

        >>> # demonstrating a failure
        >>> cmd = ShortOpts("python", m='fail')
        >>> with raises(subprocess.CalledProcessError, message="Expecting CalledProcessError"):
        ...     cmd(stderr=None)

        >>> # using append, execute python -m this
        >>> cmd = LongOpts("python").append("-m this")
        >>> cmd
        python -m this
        >>> with tmpdir.join("tmp2.txt").open("w") as fp:
        ...    cmd(stdout=fp)
        >>> with tmpdir.join("tmp2.txt").open("r") as fp:
        ...    fp.readlines()[0].strip()
        'The Zen of Python, by Tim Peters'
        """
        closeOut, closeErr = False, False
        if not stdout:
            stdout, closeOut = open(os.devnull, "w"), True
        if not stderr:
            stderr, closeErr = open(os.devnull, "w"), True

        subprocess.check_call(shlex.split(str(self)), stdout=stdout, stderr=stderr)

        if closeOut:
            stdout.close()
        if closeErr:
            stderr.close()

    def __str__(self):
        """
        :return: string representation

        >>> cmd = ShortOpts("ls", "$HOME", l="").append("-a").append("-R")
        >>> str(cmd)
        'ls $HOME -l -a -R'
        """
        return self.__repr__()

    def __repr__(self):
        """
        :return: repr

        >>> ShortOpts("ls", "$HOME", l="").append("-a").append("-R")
        ls $HOME -l -a -R
        """
        return ' '.join([str(self._exe)] +
                        [str(k) for k in self._args] +
                        [self._dash() + str(k) + (" " + str(v) if v else "")
                         for k, v in self._kwargs.items()]) + self._ext

    def _dash(self):
        return "--"


class LongOpts(CommandLine):
    def _dash(self):
        return "--"


class ShortOpts(CommandLine):
    def _dash(self):
        return "-"

