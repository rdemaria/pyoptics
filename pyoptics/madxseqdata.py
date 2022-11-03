from utils import myopen
import re

import gzip, StringIO

def myopen(fn):
  try:
    if fn.endswith('.gz'):
      return gzip.open(fn)
    else:
      return open(fn)
  except IOError:
    return StringIO.StringIO(fn)


def m(r,name,s):
  mt=r.match(s)
  return mt and [name]+list(mt.groups()) or mt

cmt=re.compile("!.*")
st=re.compile("[^;]*;")
rvar=re.compile("(^[\w\.]+)(:?=)([^;:,]+)")
ropt=re.compile("(^[\w]+)")
relem=re.compile("(^[\w\.]+):([\w\.]+)(,.+)?;")
relem2=re.compile("(^[\w\.]+)(,.+);")
relem3=re.compile("(^[\w\.]+);")

def parseattr(attr):
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

def parseseq(s):
  out=[]
  s=s.lower() #case insensitive
  s=s.replace("const ","")   #useless
  s=s.replace("shared ","")  #useless
  s=cmt.sub("",s)
  s=s.replace('\n','')
  s=s.replace(' ','')
  currseq=None
  for i in st.finditer(s):
    ss=i.group()
    tok=m(rvar,'var',ss) or m(relem,'elem',ss) \
        or m(relem2,'elem2',ss) or m(relem3,'elem3',ss)
    if not tok:
      print 'Error at %s' % ss
    elif tok[0]=='var':
      out.append(dict(type=tok[0],name=tok[1],expr=tok[3],eq=tok[2]))
    elif  tok[0] in ('elem','elem2'):
      el=dict(name=tok[1],type='elem')
      if tok[0]=='elem':
        el['parent']=tok.pop(2)
        if el['parent']=='sequence':
          currseq=el
          el['elements']=[]
      for j in parseattr(tok[2]):
        tok=m(rvar,'var',j) or m(ropt,'opt',j)
        if tok[0]=='var':
          el[tok[1]]=dict(type=tok[0],name=tok[1],expr=tok[3],eq=tok[2])
        else:
          el[tok[1]]=dict(type=tok[0],name=tok[1],expr=tok[1],eq='opt')
      if currseq:
        currseq['elements'].append(el)
      else:
        out.append(el)
    elif  tok[0] in ('elem3'):
      if tok[1]=='endsequence':
        out.append(currseq)
        currseq=None
  return out


class var:
  def __init__(self,expr,eq):
    self.expr=expr
    self.eq=eq

class elem:
  def __init__(self,**attr):
    self.__dict__.update(attr)

class madxseqdata(dict):
  def __init__(self,filename):
    s=myopen(filename).read()
    s=parseseq(s)
    for i in s:
      if i['type']=='var':
        self[i['name']]=var(i['expr'],i['eq'])
      elif i['type']=='elem':
        i.pop('type')
        self[i['name']]=elem(**i)

if __name__=='__main__':
  l=parseseq(open('ring_chrom.seq').read())
  print l
  a=madxseqdata('ring_chrom.seq')
  print a
