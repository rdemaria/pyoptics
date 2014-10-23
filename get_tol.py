import re

def read_aptol(fn):
  out={}
  for l in open(fn):
    if '!' not in l:
      ll=l.split(',')
      if len(ll)==4:
        name,v1,v2,v3=ll
        v1=v1.split(' ')[-1]
        v2=v2.split(' ')[-1]
        v3=v3.split(' ')[-2]
        out[name]=map(float,[v1,v3,v3])
  return out


def summ_tol(ttol,name):
  reg=re.compile(name)
  out=[ v for k,v in ttol.items() if reg.match(k)]
  a=[ min(i)*1e3 for i in zip(*out)]
  b=[ max(i)*1e3 for i in zip(*out)]
  c=[ sum(i)*1e3/len(out) for i in zip(*out)]
  vals=name,a[0],c[0],b[0],a[1],c[1],b[1],a[2],c[2],b[2]
  ln=['%-14s'%'element']
  ln.append( r'%s & %s & %s'%(r'r_{\rm min}',r'r_{\rm avg}',r'r_{\rm max}'))
  ln.append( r'%s & %s & %s'%(r'h_{\rm min}',r'h_{\rm avg}',r'h_{\rm max}'))
  ln.append( r'%s & %s & %s'%(r'v_{\rm min}',r'v_{\rm avg}',r'v_{\rm max}'))
  ln.append( r'//')
  print ' & '.join(ln)
  ln=['%-14s'%name]
  ln.append( r'%.3f & %.3f & %.3f'%(a[0],c[0],b[0]))
  ln.append( r'%.3f & %.3f & %.3f'%(a[1],c[1],b[1]))
  ln.append( r'%.3f & %.3f & %.3f'%(a[2],c[2],b[2]))
  ln.append( r'//')
  print ' & '.join(ln)
