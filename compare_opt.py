#!/usr/bin/python

import sys

import madlang


fn1=sys.argv[1]
fn2=sys.argv[2]


t1=madlang.open(fn1)
t2=madlang.open(fn2)

prratio=1e99;
prdiff=1e99;

if len(sys.argv)>3:
    prratio=float(sys.argv[3])
if len(sys.argv)>4:
    prdiff=float(sys.argv[4])


v1=None
v2=None

diff=[]
ratio=[]


for k in set(t1._keys()+t2._keys()):
    try:
      v2=t2._data[k]
    except KeyError:
      print "%s not found in %s"%(k,fn1)
    try:
      v1=t1._data[k]
    except KeyError:
      print "%s not found in %s"%(k,fn2)
    if type(v1) is float  and type(v2) is float:
        if abs(v1)<1e-7 or abs(v2)<1e-7:
            df=abs(v1-v2)
            if df>0:
              diff.append( ( abs(v1-v2) ,k))
        else:
            r=abs(v1/v2)
            if r<1:
                r=1/r
            if r >1:
              ratio.append( (r,k) )

diff.sort()
diff.reverse()
ratio.sort()
ratio.reverse()

if len(ratio)>0 and len(diff)>0:
  print "ratio %10.3e %s, diff %10.3e %s"%(ratio[0]+diff[0])
  if prratio <1e99:
     for v,k in ratio:
      if v> prratio:
         print "%-20s %10.3e %10.3e %10.3e"%(k,t1[k],t2[k],v)
      else:
        break


