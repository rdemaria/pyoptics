from madobj import Elem, Expr, Sequence, Namespace

from parser import parse, pyname,parses


def fromast(ast):
  root=Namespace()
  macros={}
  current_seq=None
  for st in ast:
    if st is None:
      continue
#    print st
    name,value=st
    nname=pyname(name)
    root.orig[nname]=name
    name=nname
    kind=value[0]
    if kind=='element':
      proto=value[0][1]
      if proto:
        proto=pyname(proto)
      if proto=='macro':
        out=[]
        for st in ast:
          out.append(st)
          if st==('statement', ('value', '};')):
            break
        macros[name]=out
      else:
        if name in root:
          ne=root[name]
        elif proto=='sequence':
          ne=Sequence(name)
          current_seq=ne
        else:
          #print proto
          proelem=root[proto]
          ne=Elem(proelem)
  #print macros




def load(fh):
  return fromast(parse(fh))

def loads(s):
  return fromast(parses(s))


def open(fn):
  if fn.endswith('.gz'):
    fh=gzip.open(fn)
  else:
    fh=file(fn)
  return load(fh)


