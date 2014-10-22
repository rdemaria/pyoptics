class Namespace(dict):
  def __init__(self,parent=None):
    if parent is None:
      parent=gbl
    self.parent=parent
    self.orig={}
  def __getitem__(self,k,opt=None):
    try:
      v=dict.__getitem__(self,k)
    except KeyError:
      v=self.parent[k]
    if hasattr(v,'eval'):
      v=v.eval(self)
    return v
  def getorig(self,k):
    return self[self.orig[k]]


class Elem(object):
  def __init__(self,name,keyword=None,**kwargs):
    self.name=name
    self.keyword=keyword
    if keyword is not None:
      for k,v in keyword.items():
        if k!='name':
          setattr(self,k,v)
    for k,v in kwargs.items():
        if k!='name':
          setattr(self,k,v)

class Expr(object):
  __slots__=('expr','code')
  def __init__(self,expr):
    self.expr=expr
    if expr.endswith(';'):
      self.expr=compile(expr,expr,'exec')
    else:
      self.expr=compile(expr,expr,'eval')
  def eval(self,lcl):
    return eval(self.expr,gbl,lcl)

class Sequence(object):
  def __init__(self,name,beam='default',elems=None,refer='centre'):
   if elems is None:
     self.elems=[]
   self.name=name
   self.beam=beam
   self.refer=refer


gbl={}
basenames=['hkicker', 'vkicker', 'ecollimator', 'instrument',
 'kicker', 'marker', 'monitor', 'multipole', 'octupole',
 'quadrupole', 'rbend', 'sbend', 'rcollimator', 'rfcavity',
 'sextupole', 'solenoid','tkicker','placeholder']

for e in basenames:
  gbl[e]=Elem(e)


import math
for k in dir(math):
  if not k.startswith('_'):
    gbl[k]=getattr(math,k)

