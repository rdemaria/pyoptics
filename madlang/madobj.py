from __future__ import division
from collections import namedtuple

import ast

def get_names(expr):
    tree=ast.parse(expr)
    out=[]
    toremove=[]
    for node in ast.walk(tree):
        if type(node) is ast.Attribute:
           name=_get_names(node.value).pop()
           out.append('.'.join([name,node.attr]))
           toremove.append(name)
        elif type(node) is ast.Name:
           out.append(node.id)
    for name in toremove:
       out.remove(name)
    return set(out)




class ExprUndefined(Exception):
    pass


def get_attrs(obj):
     import types
     if not hasattr(obj, '__dict__'):
         return []  # slots only
     if not isinstance(obj.__dict__, (dict, types.DictProxyType)):
         raise TypeError("%s.__dict__ is not a dictionary"
                         "" % obj.__name__)
     return obj.__dict__.keys()

def dir2(obj):
    attrs = set()
    if not hasattr(obj, '__bases__'):
        # obj is an instance
        if not hasattr(obj, '__class__'):
            # slots
            return sorted(get_attrs(obj))
        klass = obj.__class__
        attrs.update(get_attrs(klass))
    else:
        # obj is a class
        klass = obj
    for cls in klass.__bases__:
        attrs.update(get_attrs(cls))
        attrs.update(dir2(cls))
    attrs.update(get_attrs(obj))
    return list(attrs)

class Elem(object):
#  __slots__=['name','parent','_data','_ns','_orig']
  gbl={}
  def __init__(self,name=None,parent=None,**kwargs):
    self._data={}
    self.name=name
    self._ns=Elem.gbl
    self._parent=None
    if parent is not None:
        self._parent=parent.name
        self._data.update(parent._data)
    self._data.update(kwargs)
  def __repr__(self):
    def _str(k,v):
      if hasattr(v,'_get'):
          vv="%s -> %s"%(v,self[k])
      else:
          vv='%s'%(v)
      return "%-15s: %s"%(k,vv)
    if len(self._data)<60:
        data="\n%s\n"%"\n".join(_str(k,v) for k,v in sorted(self._data.items()))
    else:
        data='%d objects'%len(self._data)
    if self._parent is not None:
        data='%s,%s'%(self._parent,data)
    tmp="<%s: %s>"%(self.name,data)
    return tmp
  def __getitem__(self,k,opt=None):
    try:
      v=self._data[k]
    except KeyError:
      if self._ns is not None:
        v=self._ns[k]
      else:
        raise KeyError('%s not found'%(k,))
    if hasattr(v,'_get'):
      v=v._get(self,self.gbl)
    return v
  def __setitem__(self,k,v):
      self._data[k]=v
      if hasattr(v,'_bind'):
        v._bind(self,k)
  def __delitem__(self,k):
      v=self[k]
      del self._data[k]
  def _bind(self,obj,k):
     self._ns=obj
  def __contains__(self,k):
      return k in self._data
  def __getattribute__(self,k):
      if k.startswith('_'):
         return object.__getattribute__(self,k)
      else:
        try:
          return self[k]
        except KeyError:
          return object.__getattribute__(self,k)
#    print k, k in self.__class__.__dict__
#    if k in self.__class__.__dict__ or k.startswith('_'):
#       return object.__getattribute__(self,k)
#    elif k in self._data:
#        return self[k]
#    else:
#      raise AttributeError('%s missing in %s'%(k,self))
  def __setattr__(self,k,v):
    if k.startswith('_'):
        self.__dict__[k]=v
    else:
        self._data[k]=v
  def __delattr__(self,k):
    if k.startswith('_'):
        del self.__dict__[k]
    else:
        del self._data[k]
  def __dir__(self):
    out=dir2(self)
    out.extend(self._data.keys())
    return out
  def __iter__(self):
    return self._data.__iter__()
  def __hasattr__(self,k):
    return k in self._data
  def build_dep(self):
    out={}
    for key,att in self._data.items():
        if hasattr(att,'expr'):
            for name,idx,expr in att.get_names():
                out.setdefault(name,[]).append((key,idx,expr))
        elif hasattr(att,'_data'):
            for name,lst in att.build_dep().items():
                for ll in lst:
                  out.setdefault(name,[]).append((key,)+ll)
    return out
  def print_dep(self,key,indent='',deps=None):
    if deps is None:
      deps=self.build_dep()
    for dep in deps.get(key,[]):
        expr=dep[-1]
        idx=dep[-2]
        objs='.'.join(dep[:-2])
        if idx is not None:
            objs='%s[%d]'%(objs,idx)
        print "%s%-20s = %s"%(indent,objs,expr)
        if objs in deps:
            self.print_dep(objs,indent+'  ',deps)
  def __call__(self,expr):
      return eval(expr,{},self)


class Expr(object):
  __slots__=('expr')
  def __init__(self,expr):
    self.expr=expr
    if expr.endswith(';'):
      self.expr=compile(expr,expr,'exec')
    else:
      self.expr=compile(expr,expr,'eval')
  def _get(self,lcl,gbl={}):
    try:
         return eval(self.expr,gbl,lcl)
    except NameError as e:
        print("Warning %r undefined"%(self))
        #return ExprUndefined(e.message)
        return 0
  def __repr__(self):
      return "Expr(%r)"%(self.expr.co_filename)
  def get_names(self,ix=None):
      names=self.expr.co_names
      expr=[self.expr.co_filename]*len(names)
      return zip(names,[ix]*len(names),expr)

class ExprList(object):
  __slots__=('expr')
  def __init__(self,*expr):
    self.expr=[Expr(ex) for ex in expr]
  def _get(self,lcl,gbl={}):
    return [ ex._get(lcl,gbl) for ex in self.expr]
  def __repr__(self):
    exp=','.join(['%r'%ex.expr.co_filename for ex in self.expr])
    return 'ExprList(%s)' % exp
  def get_names(self):
      out=[]
      for ix,ex in enumerate(self.expr):
          out.extend(ex.get_names(ix))
      return out

class SeqElem(namedtuple('SeqElem','at From mech_sep slot_id')):
  @classmethod
  def from_dict(cls,data):
      return cls(*[data[k] for k in cls._fields])


class Sequence(Elem):
  _fields='at From mech_sep slot_id'.split()
  classes=dict(
    drift=namedtuple('drift',['l']),
    mult =namedtuple('mult','knl ksl hxl hyl l rel'),
    cav  =namedtuple('cav','vn f lag scav'),
    align=namedtuple('align','x y tilt'),
    block=namedtuple('block','elems'),
  )
  def __init__(self,name=None,parent=None,elems=None,**kwargs):
    Elem.__init__(self)
    if elems is None:
      self.elems=[]
  def append(self,name,elem):
    ne=Elem(name=name,
              at=elem._data.get('at'),
              From=elem._data.get('From'),
              slot_id=elem._data.get('slot_id'),
              mech_sep=elem._data.get('mech_sep'))
    self.elems.append(ne)
  def expand_struct(self,convert=classes):
    out=[]
    rest=[]
    count={}
    drift=convert['drift']
    mult =convert['mult']
    cav  =convert['cav']
    align=convert['align']
    block=convert['block']
    pos={}
    lasts=0
    drifts={}
    for el in self.elems:
        if el.From is None:
          pos[el.name]=el.at
    for el in self.elems:
      elem=self[el.name]
      s=el.at
      if el.From is not None:
          s+=pos[el.From]
      l=elem.l
      ldrift=s-l/2-lasts
      if lasts+l/2<s:
          dr=drifts.get(ldrift)
          if dr is None:
              drname='drift_%d'%len(drifts)
              dr=drift(l=ldrift)
              drifts[l]=dr
          out.append((drname,dr))
      if elem.keyword=='multipole':
          ne=mult(knl=elem.knl, ksl=elem.ksl,
                  l=elem.lrad, hxl=elem.knl[0],hyl=elem.ksl[0],rel=0)
          out.append((elem.name,ne))
      elif elem.keyword in ['marker','hmonitor','vmonitor','instrument',
                            'monitor','rcollimator']:
          ne=drift(l=0)
          out.append((elem.name,ne))
      elif elem.keyword in ['hkicker']:
          ne=mult(knl=[-elem.kick],ksl=[],
                       l=elem.lrad,hxl=elem.kick,hyl=0,rel=0)
          out.append((elem.name,ne))
      elif elem.keyword in ['vkicker']:
          ne=mult(knl=[],ksl=[elem.kick],
                       l=elem.lrad,hxl=0,hyl=elem.kick,rel=0)
          out.append((elem.name,ne))
      elif elem.keyword in ['rfcavity']:
          nvolt=elem.volt/self['beam'].pc/1000
          ne=cav(vn=nvolt,f=elem.freq*1e6,lag=elem.lag*360,scav=-1)
          out.append((elem.name,ne))
      else:
          rest.append((elem.name,elem))
      lasts=s+l/2
    return out,rest

elements=['hkicker', 'vkicker', 'ecollimator', 'instrument',
 'kicker', 'marker', 'monitor', 'multipole', 'octupole',
 'quadrupole', 'rbend', 'sbend', 'rcollimator', 'rfcavity',
 'sextupole', 'solenoid','tkicker','placeholder','drift',
 'hmonitor','vmonitor','sequence']

commands=['beam','twiss']

for e in elements:
  Elem.gbl[e]=Elem(e,keyword=e,l=0)

Elem.gbl['multipole']._data.update(knl=[0],ksl=[0])
Elem.gbl['rfcavity']._data.update(lag=0)
Elem.gbl['hkicker']._data.update(lrad=0)
Elem.gbl['vkicker']._data.update(lrad=0)
Elem.gbl['kicker']._data.update(lrad=0)

for e in commands:
  Elem.gbl[e]=Elem(e)

import math
for k in dir(math):
  if not k.startswith('_'):
    Elem.gbl[k]=getattr(math,k)

Elem.gbl['twopi']=2*math.pi

names=['proton','true','false','electron','ion']

for nn in names:
    Elem.gbl[nn]=nn

