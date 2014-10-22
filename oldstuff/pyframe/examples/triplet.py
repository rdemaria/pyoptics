from abpyutils import *
from scipy.optimize import fsolve

s=trtable('d0 q1 d1 q2 d2 q3 d3 q4 d4 end')

s.d0.l, s.d1.l, s.d2.l, s.d3.l, s.d4.l = 23, 2.7, 1, 2.9, 10
s.q1.l, s.q2.l = 9.2, 7.8
s.q1.knl1 = 122./23350*s.q1.l

def cons():
  s.q3.l, s.q4.l = s.q2.l, s.q1.l
  s.q2.knl1 = -s.q1.knl1/s.q1.l*s.q2.l
  s.q3.knl1, s.q4.knl1 = s.q2.knl1, s.q1.knl1

s.d0.betx,s.d0.bety=0.25,0.25

def myff(x,plot=True):
  (s.q1.l, s.q2.l, s.q1.knl1)=x.tolist()
  cons()
  s.track()
  return [s.q2.betymax()-s.q4.betxmax(),s.d4.alfx,s.d4.alfy]

x=array([s.q1.l,s.q2.l,s.q1.knl1])
myff(x)
x=fsolve(myff,x)
myff(x)
s.plotbeta()
