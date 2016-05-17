import unittest
import logging

class TestUtil(unittest.TestCase):

  def logSetup(self):
    self.logstream = StringIO.StringIO("")
    self.loghdlr = logging.StreamHandler(self.logstream)
    self.loghdlr.setLevel(logging.DEBUG)
    rootlog = logging.getLogger()
    try:
      self.stdhdlr = rootlog.handlers[0]
    except:
      self.stdhdlr = None
    rootlog.addHandler(self.loghdlr)
    self.oldLogLev = rootlog.getEffectiveLevel()
    rootlog.setLevel(logging.DEBUG)
    if self.stdhdlr != None:
      rootlog.removeHandler(self.stdhdlr)

  def logReset(self):
    rootlog = logging.getLogger()
    if self.stdhdlr != None:
      rootlog.addHandler(self.stdhdlr)
    rootlog.removeHandler(self.loghdlr)
    rootlog.setLevel(self.oldLogLev)
