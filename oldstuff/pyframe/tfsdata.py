import os,gzip
from numpy import array
from utils import pythonname

class tfsdata(object):
  def __init__(self,filename):
    self.data={}
    self.param={}
    if hasattr(filename,'__iter__'):
      file=filename
    else:
      if filename[-3:]=='.gz':
        file=gzip.open(filename)
      else:
        try:
          file=open(filename)
        except IOError:
          file=gzip.open(filename+'.gz')
    datalines=0
    for line in file:
      if line.strip():
        f=line.split()
        if (f[0] == '@'):  # descriptor lines
          try:
            self.param[f[1]]=conv(f[2],f[3])
          except:
            print "bad descriptor"," ".join(f)
        elif ( f[0] == '*'): # self.labels lines
          f.pop(0)
          f=[pythonname(l) for l in f]
          self.labels=f
          for l in self.labels: self.data[l]=[]
        elif (f[0] == '$'):  # type lines
          f.pop(0) ; self.types=f
        elif (f[0].startswith('#')):  # comment lines
          pass
        else :   # data lines
          f=map(conv,self.types,f)
          datalines+=1
          for l in self.labels:
            d=f.pop(0)
            self.data[l].append(d)
    for l in self.labels:
      self.data[l]=array(self.data[l])
    self.len=datalines

  def fromstring(self,string):
    self.fromfile(string.split('\n'))
    return self

def conv(t,i) :
  if ('e' in t): i=float(i)
  if ( ('s' in t) and (i[0]=='"') ): i=i[1:-1]
  if ('d' in t): i=int(i)
  return i

def pythonname(string):
  string=string.replace('[','')
  string=string.replace(']','')
  string=string.replace('.','_')
  string=string.replace('$','_')
  return string.lower()


def test():
  self=tfsdata()
#  self.fromfile('../twiss.lhcb1.data.gz')
  return True




