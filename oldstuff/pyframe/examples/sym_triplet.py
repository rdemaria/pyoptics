#!/afs/cern.ch/eng/sl/lintrack/Python-2.5/bin/python
from pyviews import *
from scipy.optimize import fsolve,fmin,fmin_bfgs

fs=trtable('d1 q1 d2 q2a d2a q2b d3 q3 d4 q4')
fs.d1.betx,fs.d1.bety=0.25,0.25
fs.l=[23,7,2.8,7,1.6,7,2.9,7,1,0]
grad=126/23348.89927
k1,k2,k3=grad,-grad,grad

alf1=100
alf2=20
mark=65
def mkfs(x):
  fs.q2a.l,fs.q3.l,fs.d3.l,fs.d4.l=x
  fs.q1.l     = fs.q3.l
  fs.q2b.l    = fs.q2a.l
  fs.q1.kn1l  = k1*fs.q1.l
  fs.q2a.kn1l = k2*fs.q2a.l
  fs.q2b.kn1l = k2*fs.q2b.l
  fs.q3.kn1l  = k3*fs.q3.l
  fs.track()
  bet1=fs.q2b.alfx>0 and fs.q2a.maxbety() or fs.q2b.maxbety()
  bet2=fs.q3.maxbetx()
  return [bet1-bet2,
          fs.q4.alfx-alf1,
          fs.q4.alfy-alf2,
          fs.q4.s-mark ]

x=array([   7.7083374 ,    9.05266681,    2.76558391,  1])
mkfs(x)

x=fsolve(mkfs,x)
print 'bet1_%2d  =%15.0f;' % (mark,fs.q4.betx)
print 'bet2_%2d  =%15.0f;' % (mark,fs.q4.bety)
print 'alf1_%2d  =%15.3f;' % (mark,fs.q4.alfx)
print 'alf2_%2d  =%15.3f;' % (mark,fs.q4.alfy)
print 'betapeak=%15.0f;' % fs.q3.maxbetx()
print 'l.mq1   =%15.3f;' % fs.q1.l
print 'l.mq2   =%15.3f;' % fs.q2a.l
print 'l.d3    =%15.3f;' % fs.d3.l
print 'l.d4    =%15.3f;' % fs.d4.l

fs.track().plot('betx bety')

fs.split(10).track().plot('betx bety')

