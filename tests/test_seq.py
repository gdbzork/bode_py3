import unittest
from tests import TestUtil
from bode.seq import Interval
from bode.seq import Sequence,SeqType
from bode.seq.bed import Bed
from bode.seq.homerPeak import HomerPeak

class TestInterval(TestUtil):

  def test_sanity(self):
    x = Interval("chr1",10,20,name="zork")
    self.assertEquals(x.name,"zork")

  def test_defaultName(self):
    x = Interval("chr1",10,20)
    self.assertEquals(x.name,"chr1:10-20")

  def test_saneInterval(self):
    x = Interval("chr1",10,20)
    self.assertEquals(x.saneInterval(),True)
    x = Interval("chr1",20,10)
    self.assertEquals(x.saneInterval(),False)
    x = Interval("chr1",20.3,10)
    self.assertEquals(x.saneInterval(),False)
    x = Interval("chr1",-20,10)
    self.assertEquals(x.saneInterval(),False)

  def test_name(self):
    x = Interval("chr1",10,20)
    x.name = "george"
    self.assertEquals(x.name,"george")

  def test_chrom(self):
    x = Interval("chr1",10,20)
    self.assertEquals(x.chrom,"chr1")
    x.chrom = "george"
    self.assertEquals(x.chrom,"george")

  def test_left(self):
    x = Interval("chr1",10,20)
    self.assertEquals(x.left,10)
    x.left = 11
    self.assertEquals(x.left,11)

  def test_right(self):
    x = Interval("chr1",10,20)
    self.assertEquals(x.right,20)
    x.right = 90
    self.assertEquals(x.right,90)

  def test_strand(self):
    x = Interval("chr1",10,20,strand="+")
    self.assertEquals(x.strand,"+")
    x.strand = "-"
    self.assertEquals(x.strand,"-")
    x = Interval("chr1",10,20)
    self.assertEquals(x.strand,".")

  def test_saneComparisons(self):
    x = Interval("chr1",20,30)
    y = Interval("chr9",10,20)
    self.assertEquals(x<y,True)
    self.assertEquals(y<x,False)
    self.assertEquals(y==x,False)
    y.chrom = "chr1"
    self.assertEquals(x<y,False)
    self.assertEquals(y<x,True)
    self.assertEquals(y==x,False)
    y.left = 20
    self.assertEquals(x<y,False)
    self.assertEquals(y<x,True)
    self.assertEquals(y==x,False)
    y.right = 30
    self.assertEquals(x<y,False)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,True)
    self.assertEquals(y==x,True)
    y.chrom = "chr1_Zork"
    self.assertEquals(x<y,True)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,False)
    self.assertEquals(y==x,False)

  def test_strandComparisons(self):
    x = Interval("chr1",20,30,strand=".")
    y = Interval("chr1",20,30,strand=".")
    self.assertEquals(x<y,False)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,True)
    self.assertEquals(y==x,True)
    y.strand = '+'
    self.assertEquals(x<y,True)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,False)
    self.assertEquals(y==x,False)
    y.strand = '-'
    self.assertEquals(x<y,True)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,False)
    self.assertEquals(y==x,False)
    x.strand = '+'
    self.assertEquals(x<y,True)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,False)
    self.assertEquals(y==x,False)
    x.strand = '-'
    self.assertEquals(x<y,False)
    self.assertEquals(y<x,False)
    self.assertEquals(x<=y,True)
    self.assertEquals(y<=x,True)
    self.assertEquals(y==x,True)

  def test_intervalRepr(self):
    x = Interval("chr1",10,20)
    self.assertEquals("%s"%x,"chr1:10-20")
    self.assertEquals(repr(x),"Interval(chr1,10,20,strand='.',name=\"chr1:10-20\")")
    x.strand = '+'
    self.assertEquals("%s"%x,"chr1:10-20+")
    self.assertEquals(repr(x),"Interval(chr1,10,20,strand='+',name=\"chr1:10-20\")")

class TestBed(TestUtil):

  def test_basicBed(self):
    x = Bed("chr1",10,20)
    self.assertEquals(x.chrom,"chr1")
    self.assertEquals(x.left,10)
    self.assertEquals(x.right,20)
    self.assertEquals(x.name,"chr1:10-20")
    self.assertEquals(x.strand,'.')
    self.assertEquals(x.score,0)

  def test_completeBed(self):
    x = Bed("chr1",10,20,name="zork",strand='+',score=99)
    self.assertEquals(x.chrom,"chr1")
    self.assertEquals(x.left,10)
    self.assertEquals(x.right,20)
    self.assertEquals(x.name,"zork")
    self.assertEquals(x.strand,'+')
    self.assertEquals(x.score,99)

  def test_bedRepr(self):
    x = Bed("chr1",10,20)
    self.assertEquals("%s"%x,"chr1\t10\t20\tchr1:10-20\t0\t.")
    self.assertEquals(repr(x),"Bed('chr1',10,20,name='chr1:10-20',score=0,strand='.')")

  def test_bedScore(self):
    x = Bed("chr1",10,20,score=99)
    self.assertEquals(x.score,99)
    x.score = 500
    self.assertEquals(x.score,500)

  def test_bedSaneInterval(self):
    x = Bed("chr1",10,20,score=99)
    self.assertEquals(x.saneInterval(),True)
    x.score = 1100
    self.assertEquals(x.saneInterval(),False)
    x.score = "zork"
    self.assertEquals(x.saneInterval(),False)

class TestSequence(TestUtil):

  def test_sequenceSanity(self):
    x = Sequence("ACGT")
    self.assertEquals(x.seq,"ACGT")
    self.assertEquals(x.name,None)
    self.assertEquals(x.seqType,SeqType.DNA)
    self.assertEquals(x.length,4)
