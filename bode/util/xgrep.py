import sys
import xlrd

class Xgrep(object):

  def __init__(self):
    self.fn = None
    self.fd = None

  def open(self,fn):
    self.fn = fn
    self.fd = xlrd.open_workbook(fn)

  def dumpRow(self,index,row):
    sys.stdout.write("%d"%(index,))
    for c in row:
      if c.ctype==xlrd.XL_CELL_EMPTY:
        sys.stdout.write("\t")
      if c.ctype==xlrd.XL_CELL_TEXT:
        sys.stdout.write("\t%s"%(c.value,))
      elif c.ctype==xlrd.XL_CELL_NUMBER:
        sys.stdout.write("\t%f"%(c.value,))
      elif c.ctype==xlrd.XL_CELL_DATE:
        t = xlrd.xldate_as_tuple(c.value)
        sys.stdout.write("\t%d-%d-%d@%d:%d"%(t[0],t[1],t[2],t[3],t[4]))
      elif c.ctype==xlrd.XL_CELL_BOOLEAN:
        sys.stdout.write("\t%s"%("T" if c.value==1 else "F",))
      elif c.ctype==xlrd.XL_CELL_ERROR:
        sys.stdout.write("\tERROR")
      elif c.ctype==xlrd.XL_CELL_BLANK:
        sys.stdout.write("\t")
    sys.stdout.write("\n")
    
  def search(self,tag):
    for sheet in self.fd.sheets():
      for i in range(0,sheet.nrows):
        row = sheet.row_slice(i)
        for c in row:
          if c.ctype==xlrd.XL_CELL_TEXT:
            if tag in c.value:
              self.dumpRow(i,row)
              break;
          elif c.ctype==xlrd.XL_CELL_NUMBER:
            repr = "%f" % (c.value,)
            if tag in repr:
              self.dumpRow(i,row)
              break;
