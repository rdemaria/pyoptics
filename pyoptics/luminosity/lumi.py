from matplotlib.pyplot import *

from numpy import sum,arange,exp,sqrt,pi,inf

import scipy.integrate
import scipy.optimize

clight=299792458
pmass=0.938272013e9;

def ftoint(u,alpha=1.0,Ax=1.0,Ay=1.0):
  """alpha=theta_c/2 * sqrt(emit/betxstar)
     Ax=sigma_z/betxstar
     Ay=sigma_z/betystar
  """
  dx=1+(Ax*u)**2
  dy=1+(Ay*u)**2
  n=exp(-u**2*(1+(alpha*Ax)**2/dx))
  d=sqrt(dx*dy)
  return n/d

def mkint(alpha=290e-6,Ax=1.0,Ay=1.0,debug=True):
  i=scipy.integrate.quad(ftoint,0,inf,args=(alpha,Ax,Ay))
  if debug:
    print("Integral tolerance: %e"%i[1])
  return 2/sqrt(pi)*i[0]


def luminosity(nb=2808,N=1.7e11,frev=clight/26658.8832,
               betx=0.55,bety=0.55,emit_n=3.75e-6,sigma_z=0.0755,pc=7e12,
               dsep=10,vcrab=0e6,r12=24,fcrab=400e6,debug=True):
  brho=pc/clight
  gamma=pc/pmass
  emit=emit_n/gamma
  sigma_x=sqrt(betx*emit)
  sigma_y=sqrt(bety*emit)
  thetac=dsep*sqrt(emit/betx)
  acrab=vcrab/(3e8*pc)*(4*pi*r12*fcrab)
  alpha=(thetac-acrab)/(2*sqrt(emit/betx))
  L=(nb*N**2*frev)/(4*pi*sigma_x*sigma_y)
  Ax=sigma_z/betx
  Ay=sigma_z/bety
  betw=alpha*sigma_z
  factor=mkint(alpha,Ax,Ay,debug=debug)
  LL=L*factor
  Fgeo=1/sqrt(1+((sigma_z*dsep)/(2*betx)**2))
  if debug:
    print("Head-on Luminosity [cm^-2 s^-1]: %e" % (L/1e4))
    print("Ext Crossing angle [urad]: %g" % (thetac*1e6))
    print("Crab Crossing angle [urad]: %g" % (acrab*1e6))
    print("beta_w: %g" % (alpha*sigma_z))
    print("Piwinski angle: %g" % (betw/betx))
    print("Loss factor : %g" % factor)
    print("Fgeo        : %g" % Fgeo)
    print("Luminosity  [cm^-2 s^-1] : %e" % (LL*1e-4))
  return LL


def piwibeta(sigma_z=0.0755,bety=0.075,dsep=10):
  def ftomin(x):
    return -luminosity(betx=x[0],bety=bety,sigma_z=sigma_z,
                      dsep=dsep,debug=False)
  betx,=scipy.optimize.fmin(ftomin,[3])
  return betx

def integrated_lumi(nb=2808,N=1.7e11,frev=clight/26658.8832,
                    betx=0.55,bety=0.55,emit_n=3.75e-6,sigma_z=0.0755,pc=7e12,
                    dsep=10, sigma_cross=0.110*1e-28,vcrab=0e6,debug=False):
  """cross section
  http://lhc-machine-outreach.web.cern.ch/lhc-machine-outreach/collisions.htm
  * inelastic (sin = 60 mbarn)
  * single diffractive (ssd = 12 mbarn)
  * elastic (sel  = 40 mbarn)
  """
  dt=1
  level_lumi=5e38
  lumi=luminosity(nb=nb,N=N,betx=betx,bety=bety,emit_n=emit_n,sigma_z=sigma_z,dsep=dsep)
  intlumi=0
  tt=0
  while N>0:
    while lumi>=level_lumi and alpha>5:
      alpha+=0.1
      lumi=luminosity(nb=nb,N=N,betx=betx,bety=bety,emit_n=emit_n,
                      sigma_z=sigma_z,alpha=alpha,debug=False)
      print(lumi)
    intlumi+=lumi*dt
    tt+=dt
    burnout=lumi*sigma_cross*dt
    N-=burnout
    print(tt,lumi,N,burnout)


def kfactor(betx,bety,dsepA,dsepB,vcrabA,vcrabB,emit_n=3.75e-6,sigma_z=0.0755,pc=7e12):
  La=luminosity(betx=betx,bety=bety,dsep=dsepA,vcrab=vcrabA,
                emit_n=emit_n,sigma_z=sigma_z,pc=pc)
  Lb=luminosity(betx=betx,bety=bety,dsep=dsepB,vcrab=vcrabB,
                emit_n=emit_n,sigma_z=sigma_z,pc=pc)
  print(La/Lb)


if __name__=='__main__':
  luminosity()
  luminosity(betx=0.15,bety=0.15)
  luminosity(betx=0.30,bety=0.075)
  luminosity(sigma_z=0.06,betx=0.30,bety=0.075,dsep=13)
  luminosity(sigma_z=0.06,betx=0.36,bety=0.10 ,dsep=13)
  luminosity(N=1.8e11,sigma_z=0.05,betx=0.25,bety=0.15,dsep=13,emit_n=2.5e-6)
  kfactor(betx=0.15,bety=0.15,dsepA=10,dsepB=10,vcrabA=6e6,vcrabB=-6e6)
  kfactor(betx=0.15,bety=0.15,dsepA=10,dsepB=10,vcrabA=6e6,vcrabB=0e6)
  kfactor(betx=0.15,bety=0.15,dsepA=10,dsepB=10,vcrabA=10e6,vcrabB=-10e6)
  luminosity(N=2.2e11,betx=0.15,bety=0.15,dsep=12.5,emit_n=2.5e-6)
  luminosity(N=3.5e11,nb=1404,betx=0.15,bety=0.15,dsep=11.4,emit_n=3.0e-6)
  luminosity(pc=4e12,N=1.7e11,nb=11390,betx=0.6,bety=0.6,dsep=13,emit_n=2.2e-6)
  luminosity(pc=4e12,N=1.6e11,nb=1390,betx=0.6,bety=0.6,dsep=8,emit_n=2.8e-6)
  luminosity(N=2.2e11,betx=0.15*2,bety=0.15/2,dsep=12.5*1.2,emit_n=2.5e-6)
  luminosity(N=3.5e11,nb=1404,betx=0.15*2,bety=0.15/2,dsep=11.4*1.2,emit_n=3.0e-6)
  luminosity(betx=0.30,bety=0.30,emit_n=1.54e-6,sigma_z=1e-9*clight/4,dsep=12,nb=2592,pc=6.5e12,N=1.24)


  def plot10():
    clf();i=121
    emit=2.5
    for dsep in [10,15]:
      subplot(i)
      for sigma_z,cl in zip([5,7.5,10],'rgb'):
        title(r'$N=1.2 \cdot 10^{11}\, d_{\rm sep}=%d \sigma$'%(dsep))
        #for by,ls in zip([0.075,0.10,0.15],['-','--','-.']):
        for by,ls in zip([0.075,0.10],['-','--']):
          bx=arange(0.05,1,0.01);
          lumi=[ luminosity(betx=bbx,bety=by,emit_n=emit*1e-6,sigma_z=sigma_z*1e-2,dsep=dsep,nb=2592,pc=7e12,N=1.2e11)/1e38 for bbx in bx]
          lbl=r'$\sigma_z=%g$ cm, $\beta_{\rm sep}=%g cm$'%(sigma_z,by*100)
          plot(bx*100,lumi,ls+cl,label=lbl)
          ylim(0,8)
          grid(True)
        legend()
        xlabel(r'$\beta_{\rm cross}$ [cm]')
        ylabel(r'peak luminosity $[10^{34}]$')
      i+=1
    savefig('flat12.png')

  plot10()

  def ftomin(x,by=0.15,emit=1.8,sep=12,bl=1):
    bx=x[0]
    return -luminosity(betx=bx,bety=by,emit_n=emit*1e-6,sigma_z=bl*1e-9*clight/4,dsep=sep,nb=2592,pc=6.5e12,N=1.24e11)

#emit=1.54
#out=[]
#for sep in [10,12,13]:
#  for bl in [0.9,1,1.1,1.2]:
#    for by in [0.05,0.10,0.15,0.20,0.25,0.30,0.40]:
#      bx=scipy.optimize.fmin(ftomin,[0.25],args=(by,emit,sep,bl))[0]
#      out.append((bx,by,emit,sep,bl))
#
#for bx,by,emit,sep,bl in out:
#  if by in [0.10,0.20,0.30,0.40] and sep in [10,12]:
#    print ("%5.3g "*4)%(bx,by,bl,sep)
#


