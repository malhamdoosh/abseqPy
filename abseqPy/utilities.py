from __future__ import print_function

import os
import sys
import subprocess
import shlex
import platform
import functools
import importlib

from abseqPy.config import ABSEQROOT, EXTERNAL_DEP_DIR
from abseqPy.config import MEM_GB


# temporarily overrides PATH variable with EXTERNAL_DEP_DIR/bin, IGBLASTDB and IGDATA (if they exist)
class PriorityPath:
    """
    this class is deprecated. It was initially used when abseqPy's pip automatically install
    3rd party dependencies. As that was phased out, this class no long has it use.
    In the future, if users want to provid a custom directory of binaries (the 3rd party bins),
    this class can support path overriding

    This class overrides the PATH, IGDATA and IGBLASTDB for abseqPy to work properly.
    PATH: will be updated to have all 3rd party dependencies (eg, clustalo, fastqc, leeHomMulti, etc ...)

    IGDATA: will be updated to the directory containing optional_data and internal_data (see IgBLAST setup page)

    IGBLASTDB: will be updated to the directory containing imgt_<species>_ig[hkl][vdj] fasta files as arguments
    to -germline_db_* in igblast

    """
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

        if sys.version_info[0] == 2:
            # escape paths containing '\' (windows path)
            subprocess.check_call(shlex.split(str(self).encode('string-escape')), stdout=stdout, stderr=stderr)
        else:
            subprocess.check_call(shlex.split(str(self).encode('unicode_escape').decode()),
                                  stdout=stdout, stderr=stderr)


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

        >>> ShortOpts("ls", "$HOME", l="").append("-a -R")
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


def disableFor(operatingSys):
    """
    a decorator that raises a NotImplementedError if this function was called on a non-supported operating system

    :param operatingSys: string, list, or tuple of string. Case insensitive.

            - win, Windows
            - mac, Darwin
            - lin, Linux

    :return: Not applicable.

    >>> import platform

    >>> # an example of disabled function
    >>> @disableFor(platform.system())
    ... def foo():
    ...     return "Your platform: " + platform.system()
    ...

    >>> with raises(NotImplementedError, message="Expecting NotImplementedError"):
    ...     foo()

    >>> # an example of enabled function
    >>> @disableFor([])
    ... def bar(message='default'):
    ...     return message
    ...

    >>> bar(message='hello')
    'hello'
    """

    if isinstance(operatingSys, str):
        systems = [operatingSys]
    elif isinstance(operatingSys, list) or isinstance(operatingSys, tuple):
        systems = operatingSys
    else:
        raise TypeError(str(operatingSys) + " is not a valid type. Expecting a string, list, or tuple.")

    # map short names to platform.system() names
    OPS = {
        "win": "Windows",
        "lin": "Linux",
        "mac": "Darwin",
        "windows": "Windows",
        "linux": "Linux",
        "darwin": "Darwin"
    }

    unsupported = []

    for s in systems:
        if s.lower() not in OPS:
            raise ValueError(str(operatingSys) + " is not a valid operating system name. "
                                                 "Only Win, Lin, and Mac is supported")
        unsupported.append(OPS[s.lower()])

    def _decorator(func):
        @functools.wraps(func)
        def _call(*args, **kwargs):
            if platform.system() in unsupported:
                raise NotImplementedError("Sorry, {} is not implemented for {}."
                                          .format(func.__name__, platform.system()))
            return func(*args, **kwargs)
        return _call
    return _decorator


def requires(package, fatal=False, stderr=sys.stderr):
    """
    a decorator that will NOT call the function if the requirements are not fulfilled. i.e. ALL packages specified in
    package MUST be importable. Importable is defined as a non-throwing importlib.import_module(package_name).

    if fatal is specified, and one of the packages cannot be imported, then an ImportException will be raised

    :param package: string, list, or tuple of package names.
                packages are case sensitive

    :param fatal: bool
                should a failure to import any of the package result in an exception?

    :param stderr: output stream.
                prints a message for the first un-importable package

    :return: Not applicable.

    >>> # function is called because os and sys is detected
    >>> @requires(['os', 'sys'])
    ... def foo(message):
    ...     import os, sys
    ...     return message
    ...
    >>> foo("No problemo")
    'No problemo'


    >>> # function is not called because ssys is not found (does not raise exception)
    >>> out = getfixture("tmpdir")
    >>> with out.join("tmp.txt").open("w") as fp:
    ...    @requires(['os', 'ssys'], stderr=fp)
    ...    def bar():
    ...        raise Exception("This function won't even be called")
    ...    bar()
    ...    print("No exception raised")
    ...
    No exception raised
    >>> with out.join("tmp.txt").open("r") as fp:
    ...     fp.readlines()[0].strip()
    "one of '['os', 'ssys']' cannot be found in your python path. Skipping 'bar' function call which depends on it."


    >>> # function is not called because ssys is not found (raise ImportError exception because it's fatal)
    >>> out = getfixture("tmpdir")
    >>> with out.join("tmp.txt").open("w") as fp:
    ...    @requires(['os', 'ssys'], stderr=fp, fatal=True)
    ...    def foobar():
    ...        raise Exception("This function won't even be called")
    ...    with raises(ImportError, message="Expecting ImportError"):
    ...        foobar()
    ...        print("This line is skipped because an exception will be raised")
    ...
    """
    if isinstance(package, str):
        packages = [package]
    elif isinstance(package, list) or isinstance(package, tuple):
        packages = package
    else:
        raise TypeError(str(package) + " is not a valid type. Expecting a string, list, or tuple.")

    def _decorator(func):
        def _hasPackage():
            try:
                for p in packages:
                    importlib.import_module(p)
                return True
            except Exception:           # python3 raises ModuleNotFound, while python2 raises ImportError ...
                if fatal:
                    # reraise
                    raise ImportError("{} is a required package for {} but is not found.".format(p, func.__name__))
                return False

        @functools.wraps(func)
        def _call(*args, **kwargs):
            if _hasPackage():
                return func(*args, **kwargs)
            else:
                print("one of '{}' cannot be found in your python path. "
                      "Skipping '{}' function call which depends on it."
                      .format(package, func.__name__), file=stderr)
        return _call
    return _decorator


def quote(string, quote="\""):
    """
    wraps string with quotes ("). Used mostly for paths with spaces
    :param string: original string
    :return: "string" with quotes(") wrapped around
    """
    return "{}{}{}".format(quote, string, quote)
