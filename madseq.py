import re
from math import *
from collections import defaultdict

from numpy import array,dtype

from collections import namedtuple



twopi=2*pi

comment=[re.compile(r"//.*"),
         re.compile("/\*.*?\*/", re.S),
         re.compile(r"!.*"),
         re.compile(r"real |const |shared ",re.I),
         re.compile(r"\s",re.S)]
statement=re.compile("[^;]*;")
variable=re.compile("(^[\w\.]+)(:?=)([^;:,]+)")
element=re.compile("(^[\w\.]+)(:([\w\.]+))?(,.+)?")
ropt=re.compile("(^[\w]+)")

def parse_attr(attr):
  ch=False
  out=[]
  if attr:
    for i in attr:
      if i=='{':
        ch=True
      if i=='}':
        ch=False
      if ch and (i==',') :
        out.append(' ')
      else:
        out.append(i)
  return out and (''.join(out[1:])).split(',') or out

def parse_variable(st):
  res=variable.match(st)
  if res:
    res=res.groups()
    cls = res[1]==':=' and 'expression' or 'variable'
    name=res[0]
    expr=res[2]
    if expr.startswith('{'):
      expr=expr[1:-1].split() # use trick in parseattr ',' -> ' '
    return [name,[cls,expr]]

def parse_element(st):
  res=element.match(st)
  if res:
    res=res.groups()
    name,proto,attr=res[0],res[2],res[3]
    value=[]
    res=['element',value]
    value.append(['proto',proto])
    if attr:
      for i in parse_attr(attr):
        value.append(parse_variable(i))
    return [name,res]


def parses(s):
  for r in comment:
    s=r.sub('',s)
  s=s.lower()
  current_sequence=None
  for st in statement.finditer(s):
    st=st.group()
    while 1:
      res=parse_variable(st)
      if res:
        break
      res=parse_element(st)
      if res:
        break
      break
    if res is None:
      res=(('class','statement'),('value',st))
    yield res

def parse(fh):
  return parses(fh.read())



base=view(dict(name='mad',l=0))

class madelem(view):
  def __init__(self,name='mad',proto=base,**kwargs):
    self._data={}
    self._proto=[proto]
    self._lcl=self
    self.name=name
    self._data.update(kwargs)
  def _bind(self,obj,name):
    self.name=name
    self._lcl=obj
    return self
  def __setitem__(self,k,v):
    self._data[k]=v
    if isinstance(v,expr):
      v.lcl=self._lcl
    elif hasattr(v,'_bind'):
      v._bind(self,k)

basenames=['hkicker', 'vkicker', 'ecollimator', 'instrument',
 'kicker', 'marker', 'monitor', 'multipole', 'octupole',
 'quadrupole', 'rbend', 'sbend', 'rcollimator', 'rfcavity',
 'sextupole', 'solenoid','tkicker','placeholder']
for i in basenames:
  base[i]=madelem(name=i,elemclass=i)




seqdata=namedtuple('seqdata','position element at From sequence')

class sequence(madelem):
  refer='centre'
  def _initseq(self):
    self._elements=[]
    self._element_dict=defaultdict(list)
  def _append(self, newelem):
    if 'From' in newelem:
      frm=newelem['From']
    else:
      frm=None
    pos=len(self._elements)
    data=seqdata(pos,newelem,newelem.at,frm,self)
    self._elements.append(data)
    self._element_dict[newelem.name].append(data)
  def get_data(self,name):
    return [ data for data in self._element_dict[name] ]
  def get_pos(self,name):
    return [ data.position for data in self.get_data(name)]
  def get_table(self,start=None,end=None,cols='name s l angle k1'):
    if start:
      start=self.get_data(start)[0].position
    else:
      start=0
    if end:
      end=self.get_data(end)[-1].position+1
    else:
      end=-1
    elems=self._elements[start:end]
    cols=cols.split()
    if 's' in cols:
      for i in elems:
          i.element.s=self.startpos(i)
    table=[  [ getattr(i.element,c,0.) for c in cols] for i in elems]
    table=zip(*table)
    table=map(array,table)
    res=namedtuple('result',cols)(*table)
    return res
  def get_elements(self,start=None,end=None):
    if start:
      start=self.get_data(start)[0].position
    else:
      start=0
    if end:
      end=self.get_data(end)[-1].position+1
    else:
      end=-1
    return [i.element for i in self._elements[start:end]]
  def startpos(self,data):
    if self.refer=='centre':
      position=data.at-data.element.l/2.
    else:
      raise self.refer + " Not implemented"
    if data.From:
      position+=self.startpos(data.sequence._element_dict[data.From][0])
    return position

def load(fh):
  return fromast(parse(fh))

def loads(fh):
  return fromast(parses(fh))

#_open=open
#def open(fn):
#  return load(_open(fn))


def fromast(ast,root=None,lcl=None,special=[],name='mad'):
  if root is None:
    root=madelem(name=name)
  if lcl is None:
    lcl=root
  current_seq=None
  root._originalnames={}
  for i,st in enumerate(ast):
    if st is None:
      continue
    name,value=st
    nname=pyname(name)
    root._originalnames[nname]=name
    name=nname
    kind=value[0]
    if kind=='variable':
      if name in special:
        value=evaluate(value[1],None)
      else:
        value=evaluate(value[1],lcl)
      setattr(root,name,value)
    elif kind=='expression':
      value=mkexpr(value[1])
      setattr(root,name,value)
    elif kind=='element':
      value=value[1] #[element, [[proto,.],[.,.]]] -> [[proto,.],[.,..]]]
      proto=value[0][1] #can be None
      if proto:
        proto=pyname(proto)
      if name=='endsequence':
        current_seq.n_elems=len(current_seq._elements)
        current_seq=None
        current_pos=0
      elif name=='return':
        pass
      else:
        if name in root: # element already defined
          ne=getattr(root,name)
        elif proto=='sequence': # special element sequence
          ne=sequence(name)
          ne._initseq()
          current_seq=ne
        else: # element not defined, therefore has proto
          proto=getattr(root,proto)
          ne=madelem(name,proto)
        setattr(root,name,ne)
        ne._parent=root
        fromast(value[1:],root=ne,lcl=lcl,special=['refer','From','apertype'])
        if current_seq is not None and current_seq is not ne:
          current_seq._append(ne)
  return root


def no_dots(x):
  return x.group().replace('.','_')

madname=re.compile(r'([a-z_][a-z_0-9\.]*)')
def pyname(n):
  n=n.lower()
  n=madname.sub(no_dots,n)
  n.replace('^','**')
  if n=='from':
    n='From'
  return n


def evaluate(value,lcl):
  if type(value) is str:
    try:
      value=pyname(value)
      if lcl:
        value=eval(value,globals(),lcl)
    except NameError:
      print 'Warning',value,'not evaluated'
  elif type(value) is list:
    value=[ evaluate(i,lcl) for i in value]
  return value

def mkexpr(value):
  if type(value) is str:
    value=pyname(value)
    value=expr(value)
  elif type(value) is list:
    value=mkexpr(repr(value))
  return value



