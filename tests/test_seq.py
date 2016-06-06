import unittest
from tests import TestUtil
from bode.seq import Interval
from bode.seq import Strand
from bode.seq import Sequence,SeqType
from bode.seq.bed import Bed
from bode.seq.homerPeak import HomerPeak

class TestInterval(TestUtil):

  def test_sanity(self):
    x = Interval("chr1",10,20,name="zork")
    self.assertEqual(x.name,"zork")

  def test_defaultName(self):
    x = Interval("chr1",10,20)
    self.assertEqual(x.name,"chr1:10-20")

  def test_saneInterval(self):
    x = Interval("chr1",10,20)
    self.assertEqual(x.saneInterval(),True)
    x = Interval("chr1",20,10)
    self.assertEqual(x.saneInterval(),False)
    x = Interval("chr1",20.3,10)
    self.assertEqual(x.saneInterval(),False)
    x = Interval("chr1",-20,10)
    self.assertEqual(x.saneInterval(),False)

  def test_name(self):
    x = Interval("chr1",10,20)
    x.name = "george"
    self.assertEqual(x.name,"george")

  def test_chrom(self):
    x = Interval("chr1",10,20)
    self.assertEqual(x.chrom,"chr1")
    x.chrom = "george"
    self.assertEqual(x.chrom,"george")

  def test_left(self):
    x = Interval("chr1",10,20)
    self.assertEqual(x.left,10)
    x.left = 11
    self.assertEqual(x.left,11)

  def test_right(self):
    x = Interval("chr1",10,20)
    self.assertEqual(x.right,20)
    x.right = 90
    self.assertEqual(x.right,90)

  def test_strand(self):
    x = Interval("chr1",10,20,strand="+")
    self.assertEqual(x.strand,Strand("+"))
    x.strand = "-"
    self.assertEqual(x.strand,Strand("-"))
    x = Interval("chr1",10,20)
    self.assertEqual(x.strand,Strand("."))

  def test_saneComparisons(self):
    x = Interval("chr1",20,30)
    y = Interval("chr9",10,20)
    self.assertEqual(x<y,True)
    self.assertEqual(y<x,False)
    self.assertEqual(y==x,False)
    y.chrom = "chr1"
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,True)
    self.assertEqual(y==x,False)
    y.left = 20
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,True)
    self.assertEqual(y==x,False)
    y.right = 30
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,True)
    self.assertEqual(y==x,True)
    y.chrom = "chr1_Zork"
    self.assertEqual(x<y,True)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,False)
    self.assertEqual(y==x,False)

  def test_strandComparisons(self):
    x = Interval("chr1",20,30,strand=".")
    y = Interval("chr1",20,30,strand=".")
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,True)
    self.assertEqual(y==x,True)
    y.strand = '+'
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,True)
    self.assertEqual(x<=y,False)
    self.assertEqual(y<=x,True)
    self.assertEqual(y==x,False)
    y.strand = '-'
    self.assertEqual(x<y,True)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,False)
    self.assertEqual(y==x,False)
    x.strand = '+'
    self.assertEqual(x<y,True)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,False)
    self.assertEqual(y==x,False)
    x.strand = '-'
    self.assertEqual(x<y,False)
    self.assertEqual(y<x,False)
    self.assertEqual(x<=y,True)
    self.assertEqual(y<=x,True)
    self.assertEqual(y==x,True)

  def test_intervalRepr(self):
    x = Interval("chr1",10,20)
    self.assertEqual("%s"%x,"chr1:10-20")
    self.assertEqual(repr(x),"Interval(chr1,10,20,strand='.',name=\"chr1:10-20\")")
    x.strand = '+'
    self.assertEqual("%s"%x,"chr1:10-20+")
    self.assertEqual(repr(x),"Interval(chr1,10,20,strand='+',name=\"chr1:10-20\")")

class TestBed(TestUtil):

  def test_basicBed(self):
    x = Bed("chr1",10,20)
    self.assertEqual(x.chrom,"chr1")
    self.assertEqual(x.left,10)
    self.assertEqual(x.right,20)
    self.assertEqual(x.name,"chr1:10-20")
    self.assertEqual(x.strand,Strand('.'))
    self.assertEqual(x.score,0)

  def test_completeBed(self):
    x = Bed("chr1",10,20,name="zork",strand='+',score=99)
    self.assertEqual(x.chrom,"chr1")
    self.assertEqual(x.left,10)
    self.assertEqual(x.right,20)
    self.assertEqual(x.name,"zork")
    self.assertEqual(x.strand,Strand('+'))
    self.assertEqual(x.score,99)

  def test_bedRepr(self):
    x = Bed("chr1",10,20)
    self.assertEqual("%s"%x,"chr1\t10\t20\tchr1:10-20\t0\t.")
    self.assertEqual(repr(x),"Bed('chr1',10,20,name='chr1:10-20',score=0,strand='.')")

  def test_bedScore(self):
    x = Bed("chr1",10,20,score=99)
    self.assertEqual(x.score,99)
    x.score = 500
    self.assertEqual(x.score,500)

  def test_bedSaneInterval(self):
    x = Bed("chr1",10,20,score=99)
    self.assertEqual(x.saneInterval(),True)
    x.score = 1100
    self.assertEqual(x.saneInterval(),False)
    x.score = "zork"
    self.assertEqual(x.saneInterval(),False)

class TestSequence(TestUtil):

  def test_sequenceSanity(self):
    x = Sequence("ACGT")
    self.assertEqual(x.seq,"ACGT")
    self.assertEqual(x.name,None)
    self.assertEqual(x.seqType,SeqType.DNA)
    self.assertEqual(x.length,4)
