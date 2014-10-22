from timeit import Timer

setup="""
from pyviews import trtable, tfstable, from_tfstable, match
from numpy import ones

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
mask=ones(len(yc),dtype=bool)

xn=x*1.1
ystart=mkfs(xn)
"""

t=Timer('xn,info=match(mkfun,mask,xn,bisec=30,maxsteps=100)',setup)
print sum(t.repeat(3,1))/3


