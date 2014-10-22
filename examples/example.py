import twiss
reload(twiss)

seq=[['d1',5,0,0,0,0],
     ['d2',5,0,0,0,0] ]

#s,betx...,bety,disp..
init=[0,0.5,0,0,0.25,0,0,1,1,1,1,1,1]

twiss.track(seq,init)

def ftosolve(x):
  print x
  lq=4
  ls=30
  kq1= x[0]
  kq2= x[1]
  seq=[['l1',ls  ,0,0,0       ,0],
       ['q1',lq  ,0,0,kq1*lq  ,0],
       ['q2',lq*2,0,0,kq2*lq*2,0],
       ['q3',lq  ,0,0,kq1*lq  ,0],
       ['l2',4   ,0,0,0       ,0]]
  init=[0,0.5,0,0,0.5,0,0,0.1,0.01,-0.2,-0.03,0,0]
  out=twiss.track(seq,init)
  end=out[-1]
  alfx=end[3]
  alfy=end[6]
  for l in out[-1:]:
    name,s,betx, alfx, mux , bety, alfy, muy=l[0:8]
    dx,dpx,dy,dpy=l[8:12]
    print ("%-6s "+"%8.6g "*7)%( name,s,betx,alfx,mux,bety,alfy,muy)
    #print ("%-6s "+"%8.6g "*5)%( name,s,dx,dpx,dy,dpy)
  return alfx,alfy-1


import scipy.optimize

x0=[0.005,-0.005]

ftosolve(x0)

x1=scipy.optimize.fsolve(ftosolve,x0)

scipy.optimize.broyden1(ftosolve,x1)

