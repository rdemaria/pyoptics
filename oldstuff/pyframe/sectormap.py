from numpy import fromstring


class sectormap:
  """Read a sectormap file from madX. Check indexing of t matrix"""
  def __init__(self,filename):
    self.elems=[]
    self.s={}
    self.k={}
    self.r={}
    self.t={}
    l=open(filename).readlines()
    for i in range(0,len(l),44):
      a,name=l[i].split()
      self.elems.append(name)
      self.s[name]=float(a)
      self.k[name]=fromstring(l[i+1],dtype=float,sep=' ')
      self.r[name]=fromstring(' '.join(l[i+2:i+8]),dtype=float,sep=' ').reshape(6,6).T
      self.t[name]=fromstring(' '.join(l[i+8:i+44]),dtype=float,sep=' ').reshape(6,6,6).T
