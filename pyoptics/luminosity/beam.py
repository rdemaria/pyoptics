from numpy import sum,arange,exp,sqrt,pi

import scipy.integrate
import scipy.optimize
import imp


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
  i=scipy.integrate.quad(ftoint,0,scipy.integrate.Inf,args=(alpha,Ax,Ay))
  if debug:
    print("Integral tolerance: %e"%i[1])
  return 2/sqrt(pi)*i[0]


class Expr(object):
  def __init__(self,expr):
    self.expr=expr
    self.co=compile(self.expr,'eval','eval')
  def __get__(self,obj,cls):
    return eval(self.co,globals(),obj.__dict__)


class Beam(object):
  def __getattribute__(self,k):
    from . import beam
    imp.reload(beam)
    ga=object.__getattribute__
    sa=object.__setattr__
    cls=ga(self,'__class__')
    sa(self,'__class__',getattr(beam,cls.__name__))
    return ga(self,k)
  brho=Expr('pc/clight')
  gamma=Expr('pc/pmass')
  burnout=Expr('lumi*cross_total/nb*2')
  pileup=Expr('lumi*cross_inelastic/frev/nb')
  def __init__(self,nb=2808,N=1.15e11,frev=clight/26658.8832,
                    betx=0.55,bety=0.55,emit_n=3.75e-6,sigma_z=0.0755,pc=7e12,
               alpha=5,cross_total=0.100*1e-28,cross_inelastic=0.078*1e-28):
    self.nb=nb
    self.N=N
    self.frev=frev
    self.betx=betx
    self.bety=bety
    self.emit_n=emit_n
    self.sigma_z=sigma_z
    self.pc=pc
    self.alpha=alpha
    self.cross_total=cross_total
    self.cross_inelastic=cross_inelastic

  def luminosity(self,debug=True,hourglass=True):
    self.emit=self.emit_n/self.gamma
    self.sigma_x=sqrt(self.betx*self.emit)
    self.sigma_y=sqrt(self.bety*self.emit)
    self.thetac=2*self.alpha*sqrt(self.emit/self.betx)
    self.L=(self.nb*self.N**2*self.frev)/(4*pi*self.sigma_x*self.sigma_y)
    self.Ax=self.sigma_z/self.betx
    self.Ay=self.sigma_z/self.bety
    self.betw=self.alpha*self.sigma_z
    self.piwi=self.betw/self.betx
    self.factorh=mkint(self.alpha,self.Ax,self.Ay,debug=debug)
    self.factornh=1/sqrt(1+self.piwi**2)
    if hourglass:
      self.lumi=self.L*self.factorh
    else:
      self.lumi=self.L*self.factornh
    if debug:
      print("Virtual Luminosity: %e" % self.L)
      print("Full Crossing angle [urad]: %g" % (self.thetac*1e6))
      print("beta_w: %g" % (self.alpha*self.sigma_z))
      print("Piwinski angle: %g" % (self.piwi))
      print("Loss factor hg: %g" % self.factorh)
      print("Loss factor no hg: %g" % self.factornh)
      print("Luminosity  [cm^-2 s^-1] : %e" % (self.lumi*1e-4))
      print("Pileup                   : %d" % (self.pileup))
      print()
    return self.lumi
  def __call__(self,**args):
    self.__dict__.update(args)
  def piwibeta(self,betxguess=3):
    def ftomin(x):
      self.betx=x[0]
      return -self.luminosity(debug=False)
    betx,=scipy.optimize.fmin(ftomin,[betxguess])
    return betx
  def lumi_f(self,vname):
    def f(x):
      setattr(self,vname,x)
      return self.luminosity(debug=False)
    return f
  def lumi_solve(self,target,vname):
    x0=getattr(self,vname)
    f=self.lumi_f(vname)
    def ftomin(x):
      return f(x)-target
    x=scipy.optimize.fsolve(ftomin,x0)
    print("%s=%g"%(vname,x))
    setattr(self,vname,x)
    return self
  def lumi_integrate(self,level,vname,maxtime=3600*6):
    dt=10
    tt=0
    intlumi=0
    out=[]
    while self.N>0 and tt<maxtime:
      self.alpha=5
      lumi=self.luminosity(debug=False)
      if lumi>=level:
        self.lumi_solve(level,'alpha')
        lumi=self.luminosity()
      self.N-=self.burnout*dt
      tt+=dt
      intlumi+=lumi*dt
      out.append([tt,lumi,self.N,self.alpha])
    return out


#b=Beam(N=2.5e11,bety=0.075,alpha=6.5)
#b.piwibeta()
#b.luminosity()
#b.lumi_solve(5e38,'N')
#n1=b.N
#b.alpha=13
#b.lumi_solve(5e38,'N')
#n2=b.N
#print (n2-n1)/b.burnout/3600

#b.luminosity()


#tt,lumi,N,alpha=zip(*b.lumi_integrate(5e38,'alpha'))


#out=[]
#for sz in arange(0.04,0.1,0.005):
#  bf=Beam(bety=0.075,alpha=6.5,sigma_z=sz)
#  bf.piwibeta()
#  lf=bf.luminosity()
#  lc=Beam(betx=0.15,bety=0.15,alpha=0,sigma_z=sz).luminosity()
#  lp=Beam(betx=0.15,bety=0.15,alpha=5,sigma_z=sz).luminosity()
#  out.append([sz,bf.betx,lc,lf,lp])
#
#
#sz,bf,lc,lf,lp=zip(*out)

#bf=Beam(emit_n=2.5e-6,pc=4e12).luminosity()
