"""Represent genomic intervals via the Interval class, and sequences
   via the Sequence class.
"""

import re

################################################################################

class Strand(object):
  """Represents the strand of a sequence.

     Strand may be top, represented by "+", bottom ("-") or none (".").
     For purposes of object creation, -1, 0 and 1 are accepted as bottom,
     none and top also.  For ordering purposes, top < none < bottom.
  """
  tagValue = {'+': -1, '.': 0, '-': 1}

  def __init__(self,tag):
    if tag in ('+',1):
      self._tag = '+'
    elif tag in ('.',0):
      self._tag = '.'
    elif tag in ('-',-1):
      self._tag = '-'
    else:
      raise ValueError("strand '%s' not in ('-','.','+')" % tag)

  def __eq__(self,other):
    return isinstance(other,type(self)) and self._tag == other._tag
  def __ne__(self,other):
    return not isinstance(other,type(self)) or self._tag != other._tag
  def __lt__(self,other):
    try:
      return self.tagValue[self._tag] < self.tagValue[other._tag]
    except AttributeError:
      return NotImplemented
  def __le__(self,other):
    try:
      return self.tagValue[self._tag] <= self.tagValue[other._tag]
    except AttributeError:
      return NotImplemented
  def __str__(self):
    return self._tag
  def __repr__(self):
    return "Strand(%s)" % self._tag
  def tag(self):
    return self._tag
  def isNone(self):
    return self._tag == '.'

class Interval(object):
  """Base class for representing genomic intervals.

     This class stores minimal genomic intervals, consisting of a chromosome,
     left and right endpoints, strand, and optionally a name.  The convention is
     0-based, and half-open, i.e. the left endpoint is the start of the
     interval, counting from 0, and the right endpoint is one past the end of
     the interval.  The expectation is that this class will be subclassed for
     various sorts of interval formats: bed, peaks, BAM, etc.

     The ``Interval`` class supports the usual comparison operators: ``<``,
     ``==`` and so on.
  """

  chromPat = re.compile("^(chr)?(\d+)(_\w+)?$")

  def __init__(self,chrom,left,right,strand=".",name=None):
    """Construct an interval object.

       :param chrom: The chromosome name.
       :type chrom: str
       :param left: The left end of the interval (counting from 0).
       :type left: int
       :param right: The right end of the interval (one position past the end).
       :type right: int
       :param name: The name of the interval (default "chr:left-right").
       :type name: str
    """
    self._chrom = chrom
    self._left = left
    self._right = right
    if strand == None:
      self._strand = Strand('.')
    else:
      self._strand = Strand(strand)
    if name == None:
      name = "%s:%d-%d" % (chrom,left,right)
    self._name = name

  def __str__(self):
    s = "%s:%d-%d" % (self._chrom,self._left,self._right)
    if not self._strand.isNone():
      s += "%s" % self._strand
    return s

  def __repr__(self):
    name = '"%s"' % (self._name,) if self._name != "" else "None"
    return "Interval(%s,%d,%d,strand='%s',name=%s)" % (self._chrom,self._left,self._right,self._strand.tag(),name)

  def __eq__(self,other):
    if other == None or not isinstance(other,type(self)):
      return False
    return self._chrom==other._chrom and self._left==other._left and self._right==other._right and self._strand==other._strand

  def __ne__(self,other):
    if other == None or not isinstance(other,type(self)):
      return True
    return self._chrom!=other._chrom or self._left!=other._left or self._right!=other._right or self._strand!=other._strand

  def _chromLT(self,c1,c2):
    mo1 = self.chromPat.match(c1)
    mo2 = self.chromPat.match(c2)
    if mo1 and mo2:
      i1 = int(mo1.group(2))
      i2 = int(mo2.group(2))
      if i1 == i2:
        t1 = mo1.group(3)
        t2 = mo2.group(3)
        return (t1 == None and t2 != None) or (t1 != None and t2 != None and t1 < t2)
      else:
        return i1 < i2
    elif mo1:
      return True
    elif mo2:
      return False
    else:
      return c1 < c2

  def __lt__(self,other):
    if self._chrom == other._chrom:
      return (self._left < other._left) or (self._left == other._left and self._right < other._right) or (self._left == other._left and self._right == other._right and self._strand < other._strand)
    else:
      return self._chromLT(self._chrom,other._chrom)

  def __le__(self,other):
    if self._chrom == other._chrom:
      return (self._left < other._left) or (self._left == other._left and self._right < other._right) or (self._left == other._left and self._right == other._right and self._strand <= other._strand)
    else:
      return self._chromLT(self._chrom,other._chrom)

  def _getName(self):
    return self._name
  def _setName(self,name):
    self._name = name
  name = property(_getName,_setName)
  """The name of the interval (get/set)."""

  def _getChrom(self):
    return self._chrom
  def _setChrom(self,chr):
    self._chrom = chr
  chrom = property(_getChrom,_setChrom)
  """The chromosome of the interval (get/set)."""

  def _getLeft(self):
    return self._left
  def _setLeft(self,l):
    self._left = l
  left = property(_getLeft,_setLeft)
  """The left endpoint of the interval (get/set)."""

  def _getRight(self):
    return self._right
  def _setRight(self,r):
    self._right = r
  right = property(_getRight,_setRight)
  """The right endpoint of the interval (get/set)."""

  def _getStrand(self):
    return self._strand
  def _setStrand(self,strand):
    self._strand = Strand(strand)
  strand = property(_getStrand,_setStrand)
  """The strand of the interval (get/set)."""

  def saneInterval(self):
    """Test the interval for sanity.

       The test checks that the left and right endpoints are ints, that
       left<=right, and that left>=0.  If left==right, the interval is of
       length 0.

       :rtype: bool
   """
    return isinstance(self._left,int) and isinstance(self._right,int) and self._left <= self._right and self._left >= 0 and (self._strand == None or isinstance(self._strand,Strand))

################################################################################

def enum(**enums):
  return type("Enum",(),enums)

SeqType = enum(DNA=1, RNA=2, PROTEIN=3, UNKNOWN=4)
"""Representation of different types of sequences.

   Legal values are:

   * ``seq.SeqType.DNA`` -- DNA sequence
   * ``seq.SeqType.RNA`` -- RNA sequence
   * ``seq.SeqType.PROTEIN`` -- protein sequence
"""

class Sequence(object):
  """Base class for representing sequences (strings of DNA/RNA/protein).

     The most basic sequence is just a DNA(RNA,protein) sequence, with no
     annotation, position, name, etc.  It is expected that this class will
     be subclassed for various types of sequences (fastq, fasta, BAM, etc).

     ``Sequence`` objects are considered equal if their sequences are identical,
     irrespective of their names.
  """

  legalDNA = "ACGTRYSWKMBDHVN"
  legalRNA = "ACGURYSWKMBDHVN"
  legalPROTEIN = "ACDEFGHIKLMNPQRSTVWY"

  def __init__(self,seq,name=None,type=None):
    """Create a sequence object.

    :param seq: A DNA (or RNA) sequence.
    :type seq: str
    :param name: The name of the sequence (if applicable).
    :type name: str
    :param type: The sort of sequence (``seq.SeqType.DNA``, ``seq.SeqType.RNA``, ``seq.SeqType.PROTEIN``).  Defaults to "best guess": if all characters are legal DNA (including ambiguity codes), DNA is assumed; if all RNA, RNA is assumed; if all protein, protein is assumed; otherwise unknown is assumed.
    :type type: SeqType
    """
    self._seq = seq
    self._name = name
    if type == None:
      self._type = self._guessType()
      self._typeGuessed = True
    else:
      self._type = type
      self._typeGuessed = False

  def __str__(self):
    return self._seq

  def __repr__(self):
    seq = '"%s"' % (self._seq,) if self._seq != None else "None"
    name = '"%s"' % (self._name,) if self._name != "" else "None"
    return "Sequence(seq=%s,name=%s)" % (seq,name)

  def __eq__(self,other):
    if other == None or not isinstance(other,type(self)):
      return False
    return self._seq == other._seq

  def __ne__(self,other):
    if other == None or not isinstance(other,type(self)):
      return True
    return self._seq != other._seq

  def _getName(self):
    return self._name
  def _setName(self,name):
    self._name = name
  name = property(_getName,_setName)
  """The name of the sequence (get/set)."""

  def _getSeq(self):
    return self._seq
  def _setSeq(self,newSeq):
    self._seq = newSeq
  seq = property(_getSeq,_setSeq)
  """The sequence itself (get/set)."""

  def _seqLength(self):
    return len(self._seq)
  length = property(_seqLength)
  """The length of the sequence (get).

     :rtype: int
  """

  def _getSeqType(self):
    return self._type
  def _setSeqType(self,newType):
    self._type = newType
  seqType = property(_getSeqType,_setSeqType)
  """The type of the sequence (get/set)."""

  def _allInSet(self,s1,s2):
    """True if all characters in s1 also in s2, false otherwise."""
    s1len = len(s1)
    allIn = True
    for i in range(0,s1len):
      if s1[i] not in s2:
        allIn = False
        break
    return allIn

  def _guessType(self):
    t = SeqType.UNKNOWN
    if self._allInSet(self._seq,self.legalDNA):
      t = SeqType.DNA
    elif self._allInSet(self._seq,self.legalRNA):
      t = SeqType.RNA
    elif self._allInSet(self._seq,self.legalPROTEIN):
      t = SeqType.PROTEIN
    return t
    
  def fasta(self,head=0,tail=0):
    """The sequence in Fasta format.

       :param head: The number of bp to trim from the head of the sequence.
       :type head: int
       :param tail: The number of bp to trim from the tail of the sequence.
       :type tail: int
       :rtype: str
    """
    slen = len(self._seq)
    tag = self._name if not self._name == none else "sequence"
    return ">%s\n%s\n" % (tag,self._seq[head:slen-tail])

  def saneSeq(self):
    """Test whether the object is more or less real-looking.

       This test returns ``True`` if the sequence is a string, the name is
       either ``None`` or a string, and the type (if guessed) is not
       ``SeqType.UNKNOWN`` or (if specified) the sequence matches the
       specified type (and returns ``False`` otherwise).

       :rtype: bool
    """
    sane = True
    if not isinstance(self._seq,str):
      sane = False
    elif self._name != None and not isinstance(self._name,str):
      sane = False
    elif self._typeGuessed and self._type==SeqType.UNKNOWN:
      sane = False
    else:
      if self._type == SeqType.DNA:
        if not self._allInSet(self._seq,self.legalDNA):
          sane = False
      elif self._type == SeqType.RNA:
        if not self._allInSet(self._seq,self.legalRNA):
          sane = False
      elif self._type == SeqType.PROTEIN:
        if not self._allInSet(self._seq,self.legalPROTEIN):
          sane = False
    return sane

################################################################################
