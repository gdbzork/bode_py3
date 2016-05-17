from ..seq import Interval

################################################################################

class Bed(Interval):
  """A class for representing ``BED`` format intervals.

     This class stores an interval from a ``BED``-formatted file.  It extends
     ``Interval`` by adding score information.  Further columns of
     the BED standard (from UCSC) are not supported.
  """

  def __init__(self,chrom,left,right,name=None,score=0,strand="."):
    """Create a ``Bed`` object.

       :param chrom: The chromosome name.
       :type chrom: str
       :param left: The left end of the interval (counting from 0).
       :type left: int
       :param right: The right end of the interval (one position past the end).
       :type right: int
       :param name: The name of the interval (default "chr:left-right").
       :type name: str
       :param score: The score of the interval (default 0).
       :type score: int
       :param strand: The strand of the interval (One of "+", "-", ".", default ".").
       :type strand: str
    """
 
    super(Bed,self).__init__(chrom,left,right,name=name,strand=strand)
    self._score = score

  def __str__(self):
    base = "%s\t%d\t%d" % (self._chrom,self._left,self._right)
    if self._name != None:
      base += "\t%s" % (self._name,)
      if self._score != None:
        base += "\t%g" % (self._score,)
        st = self._strand
        if st != None:
          base += "\t%s" % (st,)
        else:
          base += "\t."
    return base

  def __repr__(self):
    base = "Bed('%s',%d,%d" % (self._chrom,self._left,self._right)
    if self._name != None:
      base += ",name='%s'" % (self._name,)
      if self._score != None:
        base += ",score=%d" % (self._score,)
        st = self._strand
        if st != None:
          base += ",strand='%s'" % (st,)
    base += ")"
    return base

  def _getScore(self):
    return self._score
  def _setScore(self,score):
    self._score = score
  score = property(_getScore,_setScore)
  """The score of the interval (get/set).  Legal values are ints 0-1000."""

  def saneInterval(self):
    """Test whether the BED object is sane.

       Returns ``True`` if the interval is sane, the score is an int in the
       range 0-1000, and the strand is one of '+','-','.', and ``False``
       otherwise.

       :rtype: bool
    """
    sane = True
    if not super(Bed,self).saneInterval():
      sane = False
    elif not isinstance(self._score,int):
      sane = False
    elif self._score < 0 or self._score > 1000:
      sane = False
    return sane

################################################################################
