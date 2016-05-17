from ..seq import Interval

################################################################################

class HomerPeak(Interval):
  """Represent a peak from the HOMER peak caller.

     Extends ``Interval`` with the various columns of a HOMER peaks.txt file.
  """

  def __init__(self,chrom,left,right,name,strand,normalizedTagCount,focusRatio,findPeaksScore,totalTags,controlTags,foldChangeVsControl,pvalueVsControl,foldChangeVsLocal,pvalueVsLocal,clonalFoldChange):
    """Create a ``HomerPeak`` object.

       :param chrom: Chromosome name.
       :param left: Left end of interval.
       :param right: Right end of interval.
       :param name: Name of interval.
       :param strand: Strand of interval.
       :param normalizedTagCount: Normalized tag count.
       :param focusRatio: Focus ratio.
       :param findPeaksScore: findPeaks score.
       :param totalTags: total tags.
       :param controlTags: control tags (normalised to IP).
       :param foldChangeVsControl: Fold change vs control.
       :param pvalueVsControl: p-value vs control.
       :param foldChangeVsLocal: Fold change vs local.
       :param pvalueVsLocal: p-value vs local.
       :param clonalFoldChange: Clonal fold change.
    """
    super(HomerPeak,self).__init__(chrom,left,right,name)
    self._strand = strand
    self._normalizedTagCount = normalizedTagCount
    self._focusRatio = focusRatio
    self._findPeaksScore = findPeaksScore
    self._totalTags = totalTags
    self._controlTags = controlTags
    self._foldChangeVsControl = foldChangeVsControl
    self._pvalueVsControl = pvalueVsControl
    self._foldChangeVsLocal = foldChangeVsLocal
    self._pvalueVsLocal = pvalueVsLocal
    self._clonalFoldChange = clonalFoldChange

  def _getStrand(self):
    return self._strand
  def _setStrand(self,strand):
    self._strand = strand
  strand = property(_getStrand,_setStrand)
  """The strand of the peak (get/set).  Legal values are '+','-','.'."""

  def _getTagCount(self):
    return self._normalizedTagCount
  def _setTagCount(self,tagCount):
    self._normalizedTagCount = tagCount
  normalizedTagCount = property(_getTagCount,_setTagCount)
  """The (normalized) tag count of the peak (get/set)."""

  def _getFocusRatio(self):
    return self._focusRatio
  def _setFocusRatio(self,focusRatio):
    self._focusRatio = focusRatio
  focusRatio = property(_getFocusRatio,_setFocusRatio)
  """The focus ratio of the peak (get/set)."""

  def _getFindPeaksScore(self):
    return self._findPeaksScore
  def _setFindPeaksScore(self,findPeaksScore):
    self._findPeaksScore = findPeaksScore
  findPeaksScore = property(_getFindPeaksScore,_setFindPeaksScore)
  """The findPeaks score of the peak (get/set)."""

  def _getTotalTags(self):
    return self._totalTags
  def _setTotalTags(self,totalTags):
    self._totalTags = totalTags
  totalTags = property(_getTotalTags,_setTotalTags)
  """The total tag count of the peak (get/set)."""

  def _getControlTags(self):
    return self._controlTags
  def _setControlTags(self,controlTags):
    self._controlTags = controlTags
  controlTags = property(_getControlTags,_setControlTags)
  """The control tag count of the peak (get/set)."""

  def _getFoldChangeVsControl(self):
    return self._foldChangeVsControl
  def _setFoldChangeVsControl(self,foldChangeVsControl):
    self._foldChangeVsControl = foldChangeVsControl
  foldChangeVsControl = property(_getFoldChangeVsControl,_setFoldChangeVsControl)
  """The fold change vs control of the peak (get/set)."""

  def _getPvalueVsControl(self):
    return self._pvalueVsControl
  def _setPvalueVsControl(self,pvalueVsControl):
    self._pvalueVsControl = pvalueVsControl
  pvalueVsControl = property(_getPvalueVsControl,_setPvalueVsControl)
  """The p-value vs control of the peak (get/set)."""

  def _getFoldChangeVsLocal(self):
    return self._foldChangeVsLocal
  def _setFoldChangeVsLocal(self,foldChangeVsLocal):
    self._foldChangeVsLocal = foldChangeVsLocal
  foldChangeVsLocal = property(_getFoldChangeVsLocal,_setFoldChangeVsLocal)
  """The fold change vs local of the peak (get/set)."""

  def _getPvalueVsLocal(self):
    return self._pvalueVsLocal
  def _setPvalueVsLocal(self,pvalueVsLocal):
    self._pvalueVsLocal = pvalueVsLocal
  pvalueVsLocal = property(_getPvalueVsLocal,_setPvalueVsLocal)
  """The p-value vs local of the peak (get/set)."""

  def _getClonalFoldChange(self):
    return self._clonalFoldChange
  def _setClonalFoldChange(self,clonalFoldChange):
    self._clonalFoldChange = clonalFoldChange
  clonalFoldChange = property(_getClonalFoldChange,_setClonalFoldChange)
  """The clonal fold change of the peak (get/set)."""

  def saneInterval(self):
    sane = super(HomerPeak,self).saneInterval()
    if not sane:
      return sane
    if self._strand not in ('-','+','.'):
      sane = False
    elif not isinstance(self._normalizedTagCount,Number):
      sane = False
    elif not isinstance(self._focusRatio,Number):
      sane = False
    elif not isinstance(self._findPeaksScore,Number):
      sane = False
    elif not isinstance(self._totalTags,Number):
      sane = False
    elif not isinstance(self._controlTags,Number):
      sane = False
    elif not isinstance(self._foldChangeVsControl,Number):
      sane = False
    elif not isinstance(self._pvalueVsControl,Number):
      sane = False
    elif not isinstance(self._foldChangeVsLocal,Number):
      sane = False
    elif not isinstance(self._pvalueVsLocal,Number):
      sane = False
    elif not isinstance(self._clonalFoldChange,Number):
      sane = False
    return sane

################################################################################
