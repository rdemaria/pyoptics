from ir4 import tr,mkfs,mask,x,y,tm
from pyviews import *
from numpy import *


def newdir():
  v=0.5-rand(2)
  return v/sqrt(dot(v,v))

xn=x.copy()
xn,info=match(mkfs,mask,xn,maxsteps=20)
mustep=0.01

pl=tr.plotbeta()
txt=text(2.2,0,'%7.5f %7.5f' % (tm.end.mux, tm.end.muy))
out=[]
for i in range(200):
  v=newdir()
  print 'new dir'
  if len(out)>1000: break
  while info['success']:
    print tm.end.mux, tm.end.muy, 'success'
    xgood=xn.copy()
    out.append(  [tm.end.mux, tm.end.muy, xn.copy()])
    tm.end.mux+=mustep*v[0]
    tm.end.muy+=mustep*v[1]
    xn,info=match(mkfs,mask,xn,maxsteps=20,eps=1E-5,debug=False)
    tr=tr.track()
    txt.set_text('%7.5f %7.5f' % (tm.end.mux, tm.end.muy))
    pl.update()
  tm.end.mux, tm.end.muy,xn=out[-1]
  xn,info=match(mkfs,mask,xn,maxsteps=20,eps=1E-5,debug=False)

mux,muy,sol=zip(*out)
plot(mux,muy)

pl=tr.plotbeta()
for pmux,pmuy,x in out:
  y=mkfs(x)
  pl.update()

