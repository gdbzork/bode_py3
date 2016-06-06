
class BedFile(IntervalSet):
  """Represent a bed file."""

  def __init__(self):
    super(BedFile,self).__init__()
  
  def read(self):
    if self._fd == None:
      raise NoFileHandleError()
    line = self._fd.readline()
    flds = line.split()
    if len(flds) < 3:
      raise FileFormatError(self._fn,self._lineNum,"Need >= 3 fields in line.")
    bed = Bed(flds[0],int(flds[1]),int(flds[2])
    if len(flds) >= 4:
      bed.name = flds[3]
      if len(flds) >= 5:
        bed.score = int(flds[4])
        if len(flds) >= 6:
          bed.strand = flds[5]
    return bed

  def load(self,fn,lazy=True,fd=None):
    self._fn = fn
    if fd != None:
      self._fd = fd
    else:
      self._fd = open(fn,"r")
    self._nextLine = self._fd.readline().strip()
    self._header = []
    if self._current.split()[0] = "track"
      self._header.append(self._current)
      self._current = self._fd.readline().strip()

  def nexti(self):
    pass
