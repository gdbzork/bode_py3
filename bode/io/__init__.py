"""Represent sets of intervals and sequences, usually reading from and writing
   to files.
"""
import errno

################################################################################

class NoFileHandleError(IOError):
  """Error raised if trying to read from a file when no file handle exists."""
  def __init__(self):
    text = "No file descriptor to read from."
    super(NoFileHandleError,self).__init__(self,errno.ENOENT,text)

class FileFormatError(IOError):
  """Exception indicating unexpected format in file."""
  def __init__(self,fn,line=None,msg=None):
    lmsg = ": line %d" % line if line != None else ""
    rmsg =  ": '%s'" if msg != None else ""
    text = "Syntax error: %s%s%s" % (fn,line,msg)
    super(NoFileHandleError,self).__init__(self,errno.EUCLEAN,text)

################################################################################

class IntervalSet(object):

  def __init__(self,keep=True):
    """Initialize a new IntervalSet."""
    self._fn = None
    self._header = list()
    self._fd = None
    self._intervals = list()
    self._current = 0
    self._count = 0
    self._lineNum = -1

  def __iter__(self):
    """Returns iterator over elements of the set."""
    self._current = 0
    return self

  def read(self):
    """Read one object from _fd and append it to _intervals."""
    raise NotImplementedError

  def next(self):
    """Return next element of the interval set."""
    count = 3
    if self._current < len(self._intervals):
      rv = self._intervals[self._current]
      self._current += 1
    elif self._current == self._count:
      if self._fd != None and not self._fd.eof():
        pass
      return 
    
  def load(self,fn,lazy=True,fd=None):
    """Loads intervals from the given file (or file descriptor)."""
    raise NotImplementedError

  def append(self,interval):
    """Append an interval to the set."""
    raise NotImplementedError

  def save(self,fn):
    raise NotImplementedError

  def addHeader(self,line):
    self._header.append(line)
