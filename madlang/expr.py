from math import *
twopi=2*pi

class Empty(object):
  __slots__=()
  def __getitem__(self,k):
    raise KeyError, '"%s" key not found eventually' % k

Empty=Empty()

class expr(object):
  __slots__=('expr','lcl','gbl','unit','desc','name')
  def __init__(self,expr,unit=None,desc=None,lcl=Empty,gbl=globals()):
    if expr.endswith(';'):
      self.expr=compile(expr,expr,'exec')
    else:
      self.expr=compile(expr,expr,'eval')
    self.desc=desc
    self.unit=unit
    self.gbl=gbl
    self.lcl=lcl
  def _bind(self,obj,name):
    self.name=name
    self.lcl=obj
    return self
  def _get(self,obj=None,k=None):
    try:
      return eval(self.expr,self.gbl,self.lcl)
    except NameError, msg:
      print msg
      print "Expr: error in %s"%(self.expr.co_filename)
      return 0
  def __repr__(self):
    typ=self.__class__.__name__
    exp=repr(self.expr.co_filename)
    return '%s(%s)' % (typ,exp)
  def __str__(self):
    typ=self.__class__.__name__
    exp=repr(self.expr.co_filename)
    lcl=str(self.lcl)
    unit=desc=""
    if self.unit:
      unit=' %s' % repr(self.unit)
    if self.desc:
      desc=' %s' % repr(self.desc)
    return '%s%s%s' % (exp,unit,desc)



from view import view

def expr_list(data,proto):
  self=view()
  self._proto=[proto]
  for i in data.split('\n'):
    if i.strip():
      name,value,unit,desc=i.split(None,3)
      value=expr(value,unit,desc)
      setattr(self,name,value)
  return self
