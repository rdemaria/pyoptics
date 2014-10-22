from pymad import *
from matplotlib.pyplot import *

f=Frame()
f.kq1=0.005
f.kq2=-0.005
f.lq=4

f.ls=Elem(l=30)
f.q1=Elem(l=expr('lq'),  kn1l=expr('kq1*l'))
f.q2=Elem(l=expr('lq*2'),kn1l=expr('kq2*l'))
f.q3=Elem(l=expr('lq'),  kn1l=expr('kq1*l'))
f.le=Elem(l=4)

f.seq=Line('ls','q1','q2','q3','le')
f.init=Beam(betx=0.5,bety=0.5)

tt=f.seq.track(f.init)
tt.print_table('name s betx alfx mux bety alfy muy')


import scipy.optimize


def ftosolve(x):
  global nf
  f.kq1,f.kq2=x
  tt=f.seq.track(f.init)
  ftosolve.nf+=1
  return tt.le.alfx,tt.le.alfy

def ftomin(x):
  res=array(ftosolve(x))
  return sum(res**2)

def mk_dir():
  vec=random.rand(2)-0.5
  return vec/sqrt(sum(vec**2))

ftosolve.nf=0

x0=scipy.optimize.fsolve(ftosolve,array([0.005,-0.005]),xtol=1e-13)

print ftosolve.nf

# <codecell>

x0=scipy.optimize.fsolve(ftosolve,array([0.005,-0.005]),xtol=1e-13)

def mk2(lbl):
  res=[]
  for nn in range(5):
    out=[]
    dv= 0.2*10**-arange(-1,6,0.2)
    for d in dv:
      x1=x0+ mk_dir()*d
      ftosolve.nf=0
      if lbl=='ftosolve':
        x2=scipy.optimize.fsolve(ftosolve,x1,xtol=1e-13)
      else:
        x2=scipy.optimize.minimize(ftomin,x1,tol=1e-10,method=lbl).x
      print x0,x1,x2
      if sqrt(sum((x2-x0)**2))<1e-5:
        out.append(ftosolve.nf)
      else:
        out.append(nan)
    res.append(out)
  res=1.0*array(res)
  loglog(dv,mean(res,axis=0),label=lbl)
  return out

clf()
mk2('COBYLA')
mk2('Nelder-Mead')
mk2('Powell')
#mk2('Anneal')
mk2('L-BFGS-B')
mk2('TNC')
mk2('SLSQP')
mk2('ftosolve')
xlabel('distance from initial solution')
ylabel('number of function call to reach convergence')
title('Optimizer performance')
legend(loc=0)

show()



