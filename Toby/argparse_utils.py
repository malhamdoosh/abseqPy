import argparse
import gzip
import bz2
import io
import string



class MaybeCompressedFile(object):
  def __init__(self, mode='r', bufsize=-1):
    self._mode = mode
    self._bufsize = bufsize

  def bzip2_open(self, string):
    if 'w' in self._mode:
      if os.path.splitext(string)[1].lower() == 'bz2':
        return bz2.BZ2File(string, self._mode)
    elif 'r' in self._mode:
      try:
        temp = bz2.BZ2File(string, self._mode)
        temp.read(1)
        return bz2.BZ2File(string, self._mode)
      except IOError:
        pass
    return None

  def gzip_open(self, string):
    if 'w' in self._mode:
      if os.path.splitext(string)[1].lower() == 'gz':
        return gzip.open(string, self._mode)
    elif 'r' in self._mode:
      try:
        temp = gzip.open(string, self._mode)
        temp.read(1)
        return gzip.open(string, self._mode)
      except IOError:
        pass
    return None

  def __call__(self, string):
    f = None
    if string == '-':
      if 'r' in self._mode:
        f = sys.stdin
      elif 'w' in self._mode:
        f = sys.stdout
      else:
        raise ValueError('argument "-" with mode %r' % self._mode)
      if not hasattr(f, 'seekable'):
        f = io.open(f.fileno(), self._mode)
    else:
      try:
        return self.gzip_open(string) or self.bzip2_open(string) or io.open(string, self._mode, self._bufsize)
      except IOError as e:
        raise ArgumentTypeError('can\'t open "%s": %s' % (string, e))



__all__ = [
    'MaybeCompressedFile'
]
