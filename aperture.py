import os


import numpy as np
import matplotlib.pyplot as pl
from numpy import *

from pyoptics import *

#old functions
def plot_ap(apfn="temp/ap_ir5b1.tfs",nlim=30,ref=12):
  tfn=apfn.replace('ap_ir','twiss_ir')
  t=optics.open(tfn)
  ap=optics.open(apfn)
  t.ss=ap.s
  t.n1=ap.n1
  p=t.plot(x='ss',yl='n1')
  p.figure.gca().set_ylim(0,nlim)
  pl.plot(t.ss,t.ss*0+ref)
  p.figure.canvas.draw()
  return t,ap

class BeamEnvelope(object):
  def __init__(self,twiss,survey,ap,apfiles=None,offset=None):
    self.twiss=twiss
    self.survey=survey
    self.ap=ap
    self.apfiles=apfiles
    self.offset=offset
    self.energy=ap.param['energy']
    self.gamma=ap.param['gamma']
    self.emit_n=ap.param['exn']
    self.bbeat=ap.param['beta_beating']
    self.deltap=ap.param['dp_bucket_size']
    self.co=ap.param['co_radius']
    self.d_arc=ap.param['dqf']*ap.param['paras_dx']
    self.b_arc=ap.param['betaqfx']
    for k in self.apfiles:
        idx=self.twiss.apertype==k
        for name in self.twiss.name[idx]:
            idx2=self.ap.name==name
            self.ap.aper_1[idx2]=self.twiss.aper_1[idx][0]
            self.ap.aper_2[idx2]=self.twiss.aper_2[idx][0]
            self.ap.aper_3[idx2]=self.twiss.aper_3[idx][0]
            self.ap.aper_4[idx2]=self.twiss.aper_4[idx][0]
  def plot_aper_sx(self,st='k'):
      ap=self.ap
      idx=(ap.aper_1>0)&(ap.aper_1<1)
      lim=ap.aper_1[idx]
      lim2=ap.aper_1[idx]-ap.rtol[idx]-ap.xtol[idx]
      pl.plot(ap.s[idx],lim,st)
      pl.plot(ap.s[idx],lim2,st)
      pl.plot(ap.s[idx],-lim,st)
      pl.plot(ap.s[idx],-lim2,st)
      pl.ylabel('x [m]')
      pl.xlabel('s [m]')
  def plot_aper_sy(self,st='k'):
      ap=self.ap
      idx=(ap.aper_2>0)&(ap.aper_2<1)
      lim=ap.aper_2[idx]
      lim2=ap.aper_2[idx]-ap.rtol[idx]-ap.ytol[idx]
      pl.plot(ap.s[idx],lim,st)
      pl.plot(ap.s[idx],lim2,st)
      pl.plot(ap.s[idx],-lim,st)
      pl.plot(ap.s[idx],-lim2,st)
      pl.ylabel('y [m]')
      pl.xlabel('s [m]')
  def plot_beam_sx(self,nsig=None,color='b',n1=None,n2=None):
      ap=self.ap
      sig=sqrt(ap.betx*self.emit_n/self.gamma)[n1:n2]
      if nsig==None:
          nsig=ap.param['n1min']*ap.param['halo_h']/ap.param['halo_prim']
      env=sig*nsig*self.bbeat+self.co+ap.dx*self.deltap*self.bbeat
      ss=ap.s[n1:n2]
      xx=ap.x[n1:n2]
      pl.fill_between(ss,xx+env,xx-env,color=color,alpha=0.2)
      pl.fill_between(ss,xx+sig,xx-sig,color=color,alpha=0.5)
      pl.ylabel('x [m]')
      pl.xlabel('s [m]')
  def plot_beam_sy(self,nsig=None,color='b',n1=None,n2=None):
    ap=self.ap
    sig=sqrt(ap.bety*self.emit_n/self.gamma)[n1:n2]
    if nsig==None:
       nsig=ap.param['n1min']*ap.param['halo_v']/ap.param['halo_prim']
    env=sig*nsig*self.bbeat+self.co+ap.dy*self.deltap*self.bbeat
    ss=ap.s[n1:n2]
    yy=ap.y[n1:n2]
    pl.fill_between(ss,yy+env,yy-env,color=color,alpha=0.2)
    pl.fill_between(ss,yy+sig,yy-sig,color=color,alpha=0.5)
    pl.ylabel('y [m]')
    pl.xlabel('s [m]')
  def get_surv_range(self,a=None,b=None,ref=None):
    sss=self.survey
    ttt=self.twiss
    if a is None:
        idxa=0
    else:
      idxa=where(sss//a)[0][0]
    if b is None:
       idxb=-1
    else:
       idxb=where(sss//b)[0][-1]+1
    z=sss.z[idxa:idxb];x=sss.x[idxa:idxb];t=sss.theta[idxa:idxb];
    if ref is None:
        idref=0
    else:
       idref=where(sss//ref)[0][0]
    ct=cos(t[idref]);st=sin(t[idref])
    z=z-z[idref];x=x-x[idref]; t=t-t[idref]
    zn=z*ct+x*st
    xn=-z*st+x*ct
    cox=xn+ttt.x[idxa:idxb]
    coy=ttt.y[idxa:idxb]
    emitx=self.ap.param['exn']/self.ap.param['gamma']
    emity=self.ap.param['eyn']/self.ap.param['gamma']
    bx=sqrt(ttt.betx[idxa:idxb]*emitx)
    by=sqrt(ttt.bety[idxa:idxb]*emity)
    return zn,cox,coy,bx,by
  def get_pos(self,n):
    return self.ap.x[n],self.ap.y[n]
  def get_pos_deltap(self,n,pm=1):
    x,y=self.get_pos(n)
    dx=self.ap.dx[n]
    dy=self.ap.dy[n]
    return x+dx*self.deltap*pm,y+dy*self.deltap*pm
  def get_pos_tol_spec(self,n,pm=1):
    x,y=self.get_pos_deltap(n,pm=pm)
    betx=self.ap.betx[n]
    bety=self.ap.bety[n]
    dx=self.bbeat*self.d_arc*sqrt(betx/self.b_arc)*self.deltap
    dy=self.bbeat*self.d_arc*sqrt(bety/self.b_arc)*self.deltap
    co=self.co
    rtol=self.ap.rtol[n]
    xtol=self.ap.xtol[n]
    ytol=self.ap.ytol[n]
    return x,y,xtol,ytol,rtol+co+dx,rtol+co+dy
  def get_pos_tol(self,n,pm=1):
    x,y,xtol,ytol,rtolx,rtoly=self.get_pos_tol_spec(n,pm)
    return racetrack_to_polygon(x,y,xtol,ytol,rtolx,rtoly)
  def get_pos_btol(self,n,pm=1):
    x,y=self.get_pos_deltap(n,pm=pm)
    betx=self.ap.betx[n]
    bety=self.ap.bety[n]
    dx=self.bbeat*self.d_arc*self.deltap*sqrt(betx/self.b_arc)
    dy=self.bbeat*self.d_arc*self.deltap*sqrt(bety/self.b_arc)
    co=self.co
    return racetrack_to_polygon(x,y,0,0,co+dx,co+dy)
  def get_aperture(self,n):
    apertype=self.ap.apertype[n]
    if apertype=='RECTELLIPSE':
      h=self.ap.aper_1[n]
      v=self.ap.aper_2[n]
      a=self.ap.aper_3[n]
      b=self.ap.aper_4[n]
      return rectellipse_to_polygon(0,0,h,v,a,b)
    else:
      x,y=self.apfiles[apertype]
      return x,y
  def get_ap_defined(self):
    apt=self.ap.apertype
    lst=[n for n,ap in enumerate(apt) if ap=='RECTELLIPSE' or ap in self.apfiles]
    return lst
  def get_sigma(self,n):
    betx=self.ap.betx[n]
    bety=self.ap.bety[n]
    sigx=sqrt(betx*self.emit_n/self.gamma)*self.bbeat
    sigy=sqrt(bety*self.emit_n/self.gamma)*self.bbeat
    return sigx,sigy
  def get_halo_list(self,n):
    sigx,sigy=self.get_sigma(n)
    xtp,ytp=self.get_pos_tol(n,pm=1)
    xtn,ytn=self.get_pos_tol(n,pm=-1)
    xap,yap=self.get_aperture(n)
    halo=[]
    for alf in range(0,360,45):
      tmp=[]
      sig=sqrt((sigx*cos(alf*pi/180.))**2+(sigy*sin(alf*pi/180.))**2)
      for  x0,y0 in zip(np.r_[xtp,xtn],np.r_[ytp,ytn]):
        xp,yp,t1,t2,ii=intersect_ray_polygon(x0,y0,alf,xap,yap)
        tmp.append([t1/sig,t1,x0,y0,xp,yp])
      halo.append(sorted(tmp)[0]+[alf,sig])
    return halo
  def get_halo_min(self,n):
    return sorted(self.get_halo_list(n))[0]
  def get_halo_name(self,name):
    return self.ap.n1[self.get_n_name(name)]
  def get_halo_min_name(self,name):
    return self.ap.n1[self.get_n_name(name)].min()
  def get_halo_min_all(self,n1,n2):
    out=[]
    for n in self.get_ap_defined():
      if n>=n1 and n<=n2:
        print  n,self.ap.name[n]
        out.append([ n,self.ap.name[n],self.ap.s[n]]+self.get_halo_min(n))
    return out
  def plot_aperture(self,n,color='g',lbl='aperture'):
    x,y=self.get_aperture(n)
    pl.plot(x,y,color+'-',label=lbl)
    pl.title(self.ap.name[n])
  def plot_pos_btol(self,n,color='b',lbl='beam tolerances'):
    x,y=self.get_pos_btol(n,1)
    pl.plot(x,y,color+'-')
    x,y=self.get_pos_btol(n,-1)
    pl.plot(x,y,color+'-',label=lbl)
    pl.title(self.ap.name[n])
  def plot_pos_tol(self,n,color='c',lbl='beam and mech tolerances'):
    x,y=self.get_pos_tol(n,1)
    pl.plot(x,y,color+'-')
    x,y=self.get_pos_tol(n,-1)
    pl.plot(x,y,color+'-',label=lbl)
    pl.title(self.ap.name[n])
  def plot_halo_list(self,n):
    self.plot_aperture(n)
    self.plot_pos_tol(n)
    self.plot_pos_btol(n)
    for t1s,t1,x0,y0,xp,yp,alf,sig in self.get_halo_list(n):
      p,=pl.plot([x0,xp],[y0,yp],'r-')
    pl.title(self.ap.name[n])
    return p
  def get_halo(self,n,hr,hx,hy,pm=1):
    x,y,xtol,ytol,rtol,rtol=self.get_pos_tol_spec(n,pm)
    sigx,sigy=self.get_sigma(n)
    tmp=sqrt(2)*sqrt((hr - hx)*(hr - hy))
    sh=hr-hy+tmp
    sv=hr-hx+tmp
    sr=hx+hy-hr-tmp
    h=xtol+sh*sigx
    v=ytol+sv*sigy
    a=rtol+sr*sigx
    b=rtol+sr*sigy
    return racetrack_to_polygon(x,y,h,v,a,b)
  def plot_halo(self,n,halor=13,halox=12,haloy=12,color='m',lbl='halo'):
    x,y=self.get_halo(n,halor,halox,haloy,1)
    pl.plot(x,y,color+'-')
    x,y=self.get_halo(n,halor,halox,haloy,-1)
    pl.plot(x,y,color+'-',label=lbl)
    pl.title(self.ap.name[n])
  def get_n_name(self,name):
    return np.where(self.ap//name)[0]
  def plot_ap(self,nlim=30,ref=12):
    #self.twiss.plotap(self.ap,nlim=nlim)
    t=self.twiss
    ap=self.ap
    tlt=os.path.split(t.filename)[1].split('.')[0].split('_')[1]
    t.ss=ap.s
    t.n1=ap.n1
    t.plot(x='ss',yl='n1')
    t._plot.figure.gca().set_ylim(0,nlim)
    pl.plot(t.ss,t.ss*0+ref)
    t._plot.figure.canvas.draw()
  def get_n1_name(self,name):
    idx=np.where(self.ap//name)[0]
    ap=self.ap
    return zip(ap.n1[idx],ap.name[idx],ap.betx[idx],ap.bety[idx])
  def plot_halo_name(self,name,n1=None,color='m',lbl='halo',lblap=None):
    first=True
    for n in self.get_n_name(name):
      if first:
          lbl=lbl
          first=False
          lblap=lblap
      else:
          lbl=None
          lblap=None
      if n1 is None:
        nv=self.ap.n1[n]
      else:
        nv=n1
      self.plot_halo(n,nv,nv,nv,color=color,lbl=lbl)
      self.plot_aperture(n,color='k',lbl=lblap)
      self.plot_pos_tol(n,color=color,lbl=None)
      self.plot_pos_btol(n,color=color,lbl=None)
    #pl.axes().set_aspect('equal', 'datalim')
    pl.title(name)
    tv,tt=pl.xticks()
    pl.xticks(tv,map(str,tv*1000))
    tv,tt=pl.yticks()
    pl.yticks(tv,map(str,tv*1000))
    pl.xlabel('x [mm]')
    pl.ylabel('y [mm]')
    pl.grid(True)

def arc_to_polygon(a,b,t1,t2,steps=10):
  t=np.linspace(t1,t2,steps)
  return a*cos(t),b*sin(t)

def rectellipse_to_polygon(x0,y0,h,v,a,b,steps=10):
  arg1=h/float(a)
  if arg1>1:
    arg1=1;
  arg2=v/float(b)
  if arg2>1:
    arg2=1;
  t1=arccos(arg1)
  t2=arcsin(arg2)
  if t2<t1:
    xx=np.array([h,-h,-h,h,h])
    yy=np.array([v,v,-v,-v,v])
  else:
    x,y=arc_to_polygon(a,b,t1,t2,steps=steps)
    xx=np.r_[x, -x[::-1], -x,  x[::-1], x[0]]
    yy=np.r_[y,  y[::-1], -y, -y[::-1], y[0]]
  return xx+x0,yy+y0

def racetrack_to_polygon(x0,y0,h,v,a,b,steps=10):
  x,y=arc_to_polygon(a,b,0,pi/2,steps=steps)
  xx=np.r_[x+h, -x[::-1]-h, -x-h,  x[::-1]+h, x[0]+h]
  yy=np.r_[y+v,  y[::-1]+v, -y-v, -y[::-1]-v, y[0]+v]
  return xx+x0,yy+y0

def octagon_to_polygon(x0,y0,hv,d):
  #alf=arctan(2/sqrt(2)*d/h-1)
  a=2/sqrt(2)*d-hv
  if a<0:
    return None
  xx=np.array([hv,a,-a,-hv,-hv,-a, a, hv,hv])
  yy=np.array([a,hv, hv, a,-a,-hv,-hv,-a,a])
  return xx+x0,yy+y0




def intersect_ray_segment(x0,y0,alpha,x1,y1,x2,y2):
  """return interesection point between
    P0@angle and (P1,P2)
  """
  c=cos(alpha*pi/180)
  s=sin(alpha*pi/180)
  x1-=x0;x2-=x0; y1-=y0;y2-=y0
  det=s*(x2-x1)+c*(y1-y2)
  t1=(x2*y1-x1*y2)/det
  t2=(c*y1-s*x1)/det
  return x0+c*t1,y0+s*t1,t1,t2

def intersect_ray_polygon(x0,y0,alpha,x,y):
  """return interesection point between
    P0@angle and (P1,P2)
  """
  for ii in range(len(x)-1):
    x1=x[ii]; y1=y[ii]
    x2=x[ii+1]; y2=y[ii+1]
    xp,yp,t1,t2=intersect_ray_segment(x0,y0,alpha,x1,y1,x2,y2)
    if t2>=0 and t2<=1 and t1>=0:
      return xp,yp,t1,t2,ii
  return None

def polygon_to_pmesh(x0,y0,alpha,x,y):
  """return interesection point between
    (x0,y0)@angle and (x,y)
  """
  lp=len(x)
  ii=0
  xx=[]
  yy=[]
  count=-1
  for alf in alpha:
    while count<lp:
      iii=ii%lp;       iiii=(ii+1)%lp
      x1=x[iii]; y1=y[iii]
      x2=x[iiii]; y2=y[iiii]
      xp,yp,t1,t2=intersect_ray_segment(x0,y0,alf,x1,y1,x2,y2)
      if t2>=0 and t2<=1 and t1>=0:
        xx.append(xp)
        yy.append(yp)
        if count==-1:
           count=0
        break
      ii= (ii+1)%lp
      if count>0:
        count+=1
  return xx,yy




class EnvelopeOld(object):
  """
  s1=m.madtable('survey.lhcb1.data')
  s2=m.madtable('survey.lhcb2.data')
  t1=m.madtable('twiss.lhcb1.data')
  t2=m.madtable('twiss.lhcb2.data')
  ap1=m.envelope(s1,t1)
  ap2=m.envelope(s1,t1)

  hold(False)
  figure(figsize=(6,6))
  hold(True)
  plot(ap1.co[:,2],ap1.co[:,0])
  plot(ap1.xp[:,2],ap1.xp[:,0],'g',linewidth=.1)
  plot(ap1.yp2D[:,2],ap1.yp2D[:,0],'r',linewidth=.1)
  plot(ap1.xm[:,2],ap1.xm[:,0],'g',linewidth=.1)
  plot(ap1.ym2D[:,2],ap1.ym2D[:,0],'r',linewidth=.1)
  axis([-10000,10000,-15000,5000])

  t1.select(t1.pattern("IP1"),t1.name,t1.s,ap1.xsize,ap1.ysize)
  savefig('ring.eps',dpi=600)
  """

#    kbeta=1.1,      # beta beating
#    nsigma=9.5,      # 9.5
#    nemit=3.75E-6,   #3.75E-6
#    delta=1.129E-4, # RMS energy spread
#    tol=4.6E-3,  # CO=3mm + dtol=1.6mm
#    deltamax=8E-4,   # for chromaticity measurment
#    betamaxarc=180,  # maximum beta in the arcs
#    dxmaxarc=2,  # maximum beta in the arc
  def __init__(self,s,t, kbeta=1.1, nsigma=10, nemit=3.75E-6, delta=1.129E-4, tol=4.6E-3, deltamax=8.6E-4, betamaxarc=180, dxmaxarc=2, gamma=7400):
    self.co=zeros([len(s.x),3],float)
    self.xp=zeros([len(s.x),3],float)
    self.xm=zeros([len(s.x),3],float)
    self.yp=zeros([len(s.x),3],float)
    self.ym=zeros([len(s.x),3],float)
    self.yp2D=zeros([len(s.x),3],float)
    self.ym2D=zeros([len(s.x),3],float)
    self.xsize=zeros(len(s.x),float)
    self.ysize=zeros(len(s.x),float)

    for i in xrange(len(s.x)):
      vro = array([s.x[i],s.y[i],s.z[i]])
      theta,phi,psi = s.theta[i],s.phi[i],s.psi[i]
      betx,bety,dx,dy,x,y= t.betx[i],t.bety[i],t.dx[i],t.dy[i],t.x[i],t.y[i]
      thetam=array([[cos(theta) ,           0,sin(theta)],
           [          0,           1,         0],
           [-sin(theta),           0,cos(theta)]])
      phim=  array([[          1,          0,          0],
          [          0,cos(phi)   ,   sin(phi)],
          [          0,-sin(phi)  ,   cos(phi)]])
      psim=  array([[   cos(psi),  -sin(psi),          0],
          [   sin(psi),   cos(psi),          0],
          [          0,          0,          1]])
      wm=dot(thetam,dot(phim,psim))
      ex=dot(wm,array([1,0,0]))
      ey=dot(wm,array([0,1,0]))
      self.co[i]=vro+x * ex + y * ey
      emit=nemit/gamma
      dx= dxmaxarc*sqrt(betx/betamaxarc)*0.27+abs(dx)
      dy= dxmaxarc*sqrt(bety/betamaxarc)*0.27+abs(dy)
      xsize=kbeta* (  nsigma*sqrt(betx*emit) + deltamax*dx  ) + tol
      ysize=kbeta* (  nsigma*sqrt(bety*emit) + deltamax*dy  ) + tol
      self.xp[i]=self.co[i] + xsize * ex
      self.xm[i]=self.co[i] - xsize * ex
      self.yp[i]=self.co[i] + ysize * ey
      self.ym[i]=self.co[i] - ysize * ey
      self.yp2D[i]=self.co[i] + ysize * ex
      self.ym2D[i]=self.co[i] - ysize * ex
      self.xsize[i]=xsize
      self.ysize[i]=ysize



def plotenvelopex(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
  plotenvelope("ip5",t1,ap1,"mqy_4l5_b1","mqy_4r5_b1",t2,ap2,"mqy_4l5_b2","mqy_4r5_b2")
  """
  #select
  yip5=(ap1.co[t1._row_ref[ip],0]+ap2.co[t2._row_ref[ip],0])/2
  xip5=(ap1.co[t1._row_ref[ip],2]+ap2.co[t2._row_ref[ip],2])/2
  idxs1=t1._row_ref[sele]
  idxe1=t1._row_ref[eele]
  idxs2=t2._row_ref[sele2]
  idxe2=t2._row_ref[eele2]
  # start plot
  _p.hold(True)
  _p.title("Horizontal beam envelope")
  _p.xlabel(r"$z [\rm{m}]$")
  _p.ylabel(r"$x [\rm{m}]$")
  _p.grid(True)
  # closed orbit
  x1=ap1.co[idxs1:idxe1,2]-xip5
  y1=ap1.co[idxs1:idxe1,0]-yip5
  _p.plot(y1,x1,color=[0,0,1])
  x2=ap2.co[idxs2:idxe2,2]-xip5
  y2=ap2.co[idxs2:idxe2,0]-yip5
  _p.plot(y2,x2,color=[1,0,0])
  # beam1
  x1=ap1.xp[idxs1:idxe1,2]-xip5
  y1=ap1.xp[idxs1:idxe1,0]-yip5
  x2=ap1.xm[idxs1:idxe1,2]-xip5
  y2=ap1.xm[idxs1:idxe1,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='b',alpha=0.2)
  # beam2
  x1=ap2.xp[idxs2:idxe2,2]-xip5
  y1=ap2.xp[idxs2:idxe2,0]-yip5
  x2=ap2.xm[idxs2:idxe2,2]-xip5
  y2=ap2.xm[idxs2:idxe2,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='r',alpha=0.2)


def plotenvelopey(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
  plotenvelopey("ip1",t1,ap1,"mqy_4l1_b1","mqy_4r1_b1",t2,ap2,"mqy_4l1_b2","mqy_4r1_b2")
  """
  #select
  yip5=(ap1.co[t1._row_ref[ip],0]+ap2.co[t2._row_ref[ip],0])/2
  xip5=(ap1.co[t1._row_ref[ip],1]+ap2.co[t2._row_ref[ip],1])/2
  idxs1=t1._row_ref[sele]
  idxe1=t1._row_ref[eele]
  idxs2=t2._row_ref[sele2]
  idxe2=t2._row_ref[eele2]
  # start plot
  _p.hold(True)
  _p.title("Vertical beam envelope")
  _p.xlabel(r"$z [\rm{m}]$")
  _p.ylabel(r"$y [\rm{m}]$")
  _p.grid(True)
  # closed orbit
  x1=ap1.co[idxs1:idxe1,1]-xip5
  y1=ap1.co[idxs1:idxe1,0]-yip5
  _p.plot(y1,x1,color=[0,0,1])
  x2=ap2.co[idxs2:idxe2,1]
  y2=ap2.co[idxs2:idxe2,0]-yip5
  _p.plot(y2,x2,color=[1,0,0])
  # beam1
  x1=ap1.yp[idxs1:idxe1,1]-xip5
  y1=ap1.yp[idxs1:idxe1,0]-yip5
  x2=ap1.ym[idxs1:idxe1,1]-xip5
  y2=ap1.ym[idxs1:idxe1,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='b',alpha=0.2)
  # beam2
  x1=ap2.yp[idxs2:idxe2,1]-xip5
  y1=ap2.yp[idxs2:idxe2,0]-yip5
  x2=ap2.ym[idxs2:idxe2,1]-xip5
  y2=ap2.ym[idxs2:idxe2,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='r',alpha=0.2)


def plotbeamsep(s1,t1,s2,t2,ap1,ap2,eps=5E-10):
  idx1=t1// 'bbk_.*[lr]1'
  idx2=t2// 'bbk_.*[lr]1'
# plot(-ap1.co[idx1][:,0],ap1.co[idx1][:,2],'ob')
# plot(-ap2.co[idx2][:,0],ap2.co[idx2][:,2],'or')
  x=ap1.co[idx1][:,0]-ap2.co[idx2][:,0]
  y=ap1.co[idx1][:,1]-ap2.co[idx2][:,1]
  z=ap1.co[idx1][:,2]-ap2.co[idx2][:,2]
  ds1=sqrt(x**2+y**2+z**2)
  bx=maximum(t1.betx[idx1],t2.betx[idx2])
  by=maximum(t1.bety[idx1],t2.bety[idx2])
  s=sqrt(by*eps)
  sep1=ds1/s
  idx1=t1// 'bbk_.*[lr]5'
  idx2=t2// 'bbk_.*[lr]5'
#  plot(-ap1.co[idx1][:,0],ap1.co[idx1][:,2],'ob')
#  plot(-ap2.co[idx2][:,0],ap2.co[idx2][:,2],'or')
  x=ap1.co[idx1][:,0]-ap2.co[idx2][:,0]
  y=ap1.co[idx1][:,1]-ap2.co[idx2][:,1]
  z=ap1.co[idx1][:,2]-ap2.co[idx2][:,2]
  ds5=sqrt(x**2+y**2+z**2)
  bx=maximum(t1.betx[idx1],t2.betx[idx2])
  by=maximum(t1.bety[idx1],t2.bety[idx2])
  s=sqrt(by*eps)
  sep5=ds5/s
  sb=t1.s[idx1]-average(t1.s[idx1])
  _p.plot(sb,sep1,label='ip1')
  _p.plot(sb,sep5,label='ip5')
  _p.title('Beam beam separation')
  _p.xlabel('s [m]')
  _p.ylabel(r'$\sigma$')
  _p.ylim(0,20)
  _p.legend()
  _p.grid()

