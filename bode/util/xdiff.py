import sys
import xlrd

class Xdiff(object):

  def __init__(self):
    self.fna = None
    self.fda = None
    self.fnb = None
    self.fdb = None

  def diff(self,fna,fnb):
    self.fna = fna
    self.fnb = fnb
    self.fda = xlrd.open_workbook(fna)
    self.fdb = xlrd.open_workbook(fnb)
    
    nsheets = min(self.fda.nsheets,self.fdb.nsheets)
    for i in range(0,nsheets):
      sa = self.fda.sheet_by_index(i)
      sb = self.fdb.sheet_by_index(i)
      stopAtFirst = False
      if sa.nrows != sb.nrows:
        sys.stdout.write("Size mismatch: sheet %d: %d != %d\n" % (i,sa.nrows,sb.nrows))
        sys.stdout.write("(Stopping at first mismatch.)\n")
        stopAtFirst = True
      for j in range(0,min(sa.nrows,sb.nrows)):
        ra = sa.row_slice(j)
        rb = sb.row_slice(j)
        if len(ra) != len(rb):
          sys.stdout.write("%d: length mismatch: %d != %s\n" % (j,len(ra),len(rb)))
          if stopAtFirst:
            break
        breaking = False
        for k in range(0,min(len(ra),len(rb))):
          if ra[k].value != rb[k].value:
            sys.stdout.write("%d,%d: %s != %s\n" % (j,k,ra[k].value,rb[k].value))
            if stopAtFirst:
              breaking = True
        if breaking:
          break
