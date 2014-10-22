from pyviews import *
from numpy import *

t=trtable('q1 d1 q2')



def mk(x):
  t.d1.l=3
  t.q1.l=1
  t.q2.l=1
  t.q1.kn1l=.12*t.q1.l
  t.q2.kn1l=-.12*t.q2.l
  t.d1.kn0l=0.1*t.d1.l

mk(0)
#t.track()

tf=trtable.cat(t,-t)
tf.ptrack()
p=tf.plotbeta()
p.update()

tn=trtable.cat(t,t,-t,-t)

tn.ptrack().plotbeta()

