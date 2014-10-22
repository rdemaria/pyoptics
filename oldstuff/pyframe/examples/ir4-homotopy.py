from pyviews import *
from scipy.optimize import *

tf=tfstable('twiss.ir4b1.data')
tr=from_tfstable(tf)
tr.track()

quads="mqt_13l4_b1 mqt_12l4_b1 mqtli_11l4_b1 mqml_10l4_b1 mqmc_9l4_b1 mqm_9l4_b1 mqml_8l4_b1 mqm_7l4_b1 mqy_6l4_b1 mqy_5l4_b1 mqy_5r4_b1 mqy_6r4_b1 mqm_7r4_b1 mqml_8r4_b1 mqmc_9r4_b1 mqm_9r4_b1 mqml_10r4_b1 mqtli_11r4_b1 mq_12r4_b1 mq_13r4_b1, kn1l"

cons='ip4, dx dpx; end, betx bety alfx alfy dx dpx mux muy'

idx=tr.getaddr(quads)
idy=tr.getaddr(cons)

x=tr._all[idx]

def mkfs(x):
  tr._all[idx]=x
  tr.track()
  return tr._all[idy]

def mkfun(x):
  return mkfs(x)-yc

ytarget=mkfs(x)
yc=ytarget

seed=rand(len(x))
xn=x*(1+0.4*seed)
ystart=mkfs(xn)

pl=tr.plotbeta()
lstep=0.01
l=0.
for i in range(1000):
  if abs(l-1)<1E-10:
    print 'success'
    break
  l+=lstep
  yc=ystart + l * (ytarget-ystart)
  xn,info=match(mkfun,xn,bisec=50,maxsteps=50,dbg=False)
  a=mkfs(xn)
  pl.update()
  if not info['success']:
    print 'restart'
    xn=xn*(1+3*lstep*rand(len(x)))
    ystart=mkfs(xn)
    l=0





def newdir():
  v=rand(2)
  return v/sqrt(dot(v,v))


#pl=tr.plotbeta()


