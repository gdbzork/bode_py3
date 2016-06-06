import sys
import os
import random
import subprocess
import shutil
import string

#
# Pipeline to run motif programs, from summits to results
#
DEFAULTS = {
  "genome":"/data02/genomes/Homo_sapiens_GRC_h37/hsa.fa",
  "flank":100,
  "sample":400,
  "meme_opts":["-dna","-mod","anr","-nmotifs","5","-maxsites","400","-minw","6","-maxw","16","-revcomp"],
  "siomics_opts":["-w","10","-m","5","-r","20","-s","40","-c","0.01"],
  "siomics_path":"/usr/local/src/SIOMICS_V1.3_Linux",
  "tomtom_targets":"/home/brown22/refdata/matrices/transfac_matrix_meme.dat",
  "meme-chip_opts":["-run-mast","-run-ama","-meme-minw","8","-meme-maxw","15"],
  "homer_opts":["-size","200","-mask"],
  "homer_genome":"hg19"
}
def getDefaults():
  return dict(DEFAULTS)

# Step: Sort summits, sample from top half
def sample_summits(inFN,outFN,conf):
  # read summits
  fd = open(inFN)
  d = list()
  for line in fd:
    flds = line.split()
    d.append((float(flds[4]),flds[0],int(flds[1]),flds[3]))
  fd.close()
  d.sort()
  d.reverse()
  dlen = len(d) / 2
  if dlen < conf["sample"]:
    dlen = min(len(d),conf["sample"])
  d = d[0:dlen]
  if dlen > conf["sample"]:
    random.shuffle(d)
  outFD = open(outFN,"w")
  for i in range(0,min(dlen,conf["sample"])):
    (score,chrom,pos,name) = d[i]
    outFD.write("%s\t%s\t%s\t%s\t%s\t+\n" % (chrom,pos,pos+1,name,score))
  outFD.close()

# Step: Expand reads into wider regions
def expand_reads(inFN,outFN,conf):
  inFD = open(inFN)
  outFD = open(outFN,"w")
  flank = conf["flank"]
  for line in inFD:
    flds = line.split()
    p = int(flds[1])
    l = max(p-flank,0)
    r = p+flank
    if len(flds) > 5:
      outFD.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (flds[0],l,r,flds[3],flds[4],flds[5]))
    else:
      outFD.write("%s\t%s\t%s\t%s\t%s\t+\n" % (flds[0],l,r,flds[3],flds[4]))
  inFD.close()
  outFD.close()

# Step: Generate Fasta
def bed2fasta(inFN,outFN,conf):
  cmd = ["/usr/local/bin/bedtools","getfasta","-fi",conf["genome"],"-fo",outFN,"-bed",inFN]
  proc = subprocess.Popen(cmd)
  proc.wait()

# Step: Mask Fasta (lowercase -> N)
def maskFasta(inFN,outFN):
  inFD = open(inFN)
  outFD = open(outFN,"w")
  trans = string.maketrans("acgt","NNNN")
  for line in inFD:
    if line[0] != '>':
      line = line.translate(trans)
    outFD.write(line)
  inFD.close()
  outFD.close()

# Step: Run MEME
def run_meme(inFN,outDir,conf):
  cmd = ["/usr/local/bin/meme",inFN,"-oc",outDir] + conf["meme_opts"]
  proc = subprocess.Popen(cmd)
  proc.wait()

# Step: Run MEME-CHIP
def run_memechip(inFN,outDir,conf):
  cmd = ["/usr/local/bin/meme-chip","-oc",outDir,"-db",conf["tomtom_targets"]] + conf["meme-chip_opts"] + [inFN]
  proc = subprocess.Popen(cmd)
  proc.wait()

# Step: Run HOMER
def run_homer(inFN,outDir,conf,length):
  cmd = ["/data02/jclab/homer/bin/findMotifsGenome.pl",inFN,conf["homer_genome"],outDir] + conf["homer_opts"] + ["-len",length]
  proc = subprocess.Popen(cmd)
  proc.wait()

# Step: Run TomTom
def run_tomtom(inFN,outDir,conf):
  cmd = ["/usr/local/bin/tomtom","-oc",outDir,inFN,conf["tomtom_targets"]]
  proc = subprocess.Popen(cmd)
  proc.wait()

# Step: Run SIOMICS
def run_siomics(fasta,memedir,conf):
  realdir = os.getcwd()
  fullpath = os.path.join(realdir,fasta)
  os.chdir(conf["siomics_path"])
  shutil.rmtree(memedir,ignore_errors=True)
  try:
    os.unlink(fasta)
  except OSError:
    pass
  os.symlink(fullpath,fasta)
  cmd = ["python","SIOMICS.py","-i",fasta,"-o",memedir] + conf["siomics_opts"]
  sys.stderr.write("cmd: %s\n" % (" ".join(cmd),))
  proc = subprocess.Popen(cmd)
  proc.wait()
  shutil.copytree(memedir,os.path.join(realdir,memedir))
  os.chdir(realdir)

def meme(inFN,conf):
  sampled = inFN.replace("summits","sampled")
  sample_summits(inFN,sampled,conf)
  flanks = inFN.replace("summits","flanks")
  expand_reads(sampled,flanks,conf)
  fasta = flanks.replace("bed","fa")
  bed2fasta(flanks,fasta,conf)
  memedir = fasta.replace(".fa","_MEME")
  run_meme(fasta,memedir,conf)
  tomtom = memedir.replace("MEME","TOMTOM")
  run_tomtom(os.path.join(memedir,"meme.xml"),tomtom,conf)

def siomics(inFN,conf):
  sampled = inFN.replace("summits","sampled")
  sample_summits(inFN,sampled,conf)
  flanks = inFN.replace("summits","flanks")
  expand_reads(sampled,flanks,conf)
  fasta = flanks.replace("bed","fa")
  bed2fasta(flanks,fasta,conf)
  memedir = fasta.replace(".fa","_SIOMICS")
  run_siomics(fasta,memedir,conf)

def meme_chip(inFN,conf):
  flanks = inFN.replace("summits","flanks")
  expand_reads(inFN,flanks,conf)
  fasta = flanks.replace("bed","fa")
  bed2fasta(flanks,fasta,conf)
  masked = fasta.replace("flanks","masked")
  maskFasta(fasta,masked)
  memedir = masked.replace(".fa","_MEMECHIP")
  run_memechip(masked,memedir,conf)

def homer(inFN,conf,length):
  flanks = inFN.replace("summits","flanks")
  expand_reads(inFN,flanks,conf)
  homerdir = flanks.replace(".bed","_HOMER")
  run_homer(flanks,homerdir,conf,length)
