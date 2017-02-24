import re

import matplotlib.pyplot as pl

from numpy import *

from pydataobj import dataobj
import tfsdata

from poly_fit import poly_fit, poly_print, poly_val

class StrTable(dataobj):
  scale=23348.89927
  @classmethod
  def open(cls,fn):
    obj=cls(tfsdata.open(fn))
    return obj
  def get_vars(self,reg):
    rxp=re.compile(reg)
    return sorted(l for l in self.keys() if rxp.match(l))
  def get_kq(self,n):
    return self.get_vars(r'kq[xt]?l?%da?\.'%n)
  def get_phases(self):
    out=self.get_vars(r'mu[xy]ip[1-8]b[12]$')
    out+=self.get_vars(r'mu[xy]ip[1-8]b[12]_l')
    return out
  def get_acb(self,n,knob='on_sep'):
    out=[]
    for n in self.get_vars('acb.*%d\.[lr][1-8].*%s'%(n,knob)):
      if sum(abs(self[n]))>0:
        out.append(n)
    return sorted(out)
  def plot_acb(self,n,knob,n1,n2,x=None,scale=1,brho=None):
    if brho is None:
        scale*=1e6
    else:
        scale*=brho
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    ks=self.get_acb(n,knob)
    for k in ks:
      pl.plot(xv,self[k][n1:n2]*scale,label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    if  brho is None:
      pl.ylabel(r'angle [$\mu$rad]')
    else:
      pl.ylabel('k0l [Tm]')
  def get_triplet(self,trim=True):
    tp=self.get_vars('kqx[123]?\.[rl][2815]')
    if trim:
        tp+=self.get_vars('kqtx[123]\.[rl][2815]')
    return tp
  def plot_triplet(self,n1,n2,x=None):
    scale=self.scale
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    ks=self.get_triplet()
    for k in ks:
      pl.plot(xv,abs(self[k][n1:n2]*scale),label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel('k [T/m]')
    a,b=pl.xticks()
    pl.xticks(a[::2])
  def plot_2in1(self,kq,n1,n2,x=None,sign=False,ylab='k [T/m]'):
    scale=self.scale
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    for k in self.get_kq(kq):
      kv=self[k][n1:n2]*scale
      if sign or kv[0]>0:
        pl.plot(xv,kv,label=k)
      else:
        pl.plot(xv,-kv,label='-'+k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(ylab)
    a,b=pl.xticks()
    pl.xticks(a[::2])
  def plot_ipbeta(self,n1=None,n2=None,x=None):
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    for k in self.get_vars('bet'):
      kv=self[k][n1:n2]
      pl.plot(xv,kv,label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(r"$\beta$ [m]")
  def plot_phase(self,n1=None,n2=None,x=None):
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    colors='bgrc'
    for k in self.get_phases():
      kv=self[k][n1:n2]
      pl.plot(xv,kv,color=colors[0],label=k)
      colors=colors[1:]+colors[0]
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(r'mu [2$\pi$]')
    a,b=pl.xticks()
    pl.xticks(a[::2])
  def plot_squeeze(self,n1=0,n2=None,x=None):
    fig=pl.figure('squeeze',figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    if len(self.get_vars('kqx'))>0:
      pl.subplot(3,4,1)
      self.plot_triplet(n1,n2,x=x)
    for n in range(4,11):
      if len(self.get_kq(n))>0:
        pl.subplot(3,4,n-2)
        self.plot_2in1(n,n1,n2,x=x,sign=False)
    for n in range(11,14):
      if len(self.get_kq(n))>0:
        pl.subplot(3,4,n-2)
        self.plot_2in1(n,n1,n2,x=x,sign=True)
    pl.subplot(3,4,12)
    #self.plot_ipbeta(n1,n2,x=x)
    self.plot_phase(n1,n2,x=x)
    #pl.tight_layout()
    self.xvar=x
    return self
  def plot_betsqueeze(self,n1=0,n2=None,figname=None):
    x=self.get_vars('betxip')[0]
    if figname is None:
      fig=pl.figure(x,figsize=(16,12))
    else:
      fig=pl.figure(figname,figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    pl.subplot(3,4,1)
    if len(self.get_vars('kqx'))>0:
      self.plot_triplet(n1,n2,x=x)
    for n in range(4,11):
      pl.subplot(3,4,n-2)
      self.plot_2in1(n,n1,n2,x=x,sign=False)
    for n in range(11,14):
      pl.subplot(3,4,n-2)
      self.plot_2in1(n,n1,n2,x=x,sign=True)
    pl.subplot(3,4,12)
    self.plot_phase(n1,n2,x=x)
    pl.tight_layout()
    self.xvar=x
    return self
  def plot_knobs(self,n1=0,n2=None,figname=None,scales=[1,1]):
    x=self.get_vars('betxip')[0]
    if figname is None:
      fig=pl.figure('knobs',figsize=(16,12))
    else:
      fig=pl.figure(figname,figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    for ii,(knob,scale) in enumerate(zip(['on_x','on_sep'],scales)):
       pl.subplot(2,3,1+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(1,knob,n1,n2,x=x,scale=scale)
       self.plot_acb(2,knob,n1,n2,x=x,scale=scale)
       self.plot_acb(3,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
       pl.subplot(2,3,2+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(4,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
       pl.subplot(2,3,3+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(5,knob,n1,n2,x=x,scale=scale)
       if self.get_acb(6):
         self.plot_acb(6,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
    pl.tight_layout()
    self.xvar=x
    return self
  def linint(t,name,val1,val2,x,par='ttt'):
    n1=where(t[x]==val1)[0][-1]
    n2=where(t[x]==val2)[0][-1]
    tmp= "%s:=(%g)*(1-%s)+(%g)*(%s);"
    print tmp%(name,t[name][n1],par,t[name][n2],par)
  def button_press(self,event):
    self.event=event
  def poly_fit(self,var,order,n1=None,n2=None,param=None,slope0=[]):
    if param is None:
        param=[k for k in self.keys() if k.startswith('betx')][0]
    scale=self.scale
    x=self[param][n1:n2]; y=self[var][n1:n2]
    x0=[x[0],x[-1]]
    y0=[y[0],y[-1]]
    xp0=[]; yp0=[]
    for idx in slope0:
        xp0.append(x[idx])
        yp0.append(0)
    pol=poly_fit(order,x,y,x0,y0,xp0,yp0)
    out="%s:=%s;"%(var, poly_print(pol,x=param,power='^'))
    if re.match('kt?qx[0-9]?[ab]?\.',var):
        n=3
    else:
      n=re.match('[kqtxl]+([0-9]+)\.',var)
      if n is None:
         n=14;scale=1
      else:
         n=int(n.groups()[0])
    pl.subplot(3,4,n-2)
    xv=self[param]
    yv=poly_val(pol,xv)
    print '!',' '.join(['%2d'%i for i in  sign(diff(yv))])
    if n<11:
        yv=abs(yv)
    pl.plot(xv,yv*scale)
    return out
  def poly_fit_all(self,order,param,fn,n1=None,n2=None):
    out=[]
    for kq in self.get_triplet():
        out.append(self.poly_fit(kq,order,param=param,n1=n1,n2=n2))
    for n in range(4,14):
        for kq in self.get_kq(n):
            out.append(self.poly_fit(kq,order,param=param,n1=n1,n2=n2))
    open(fn,'w').write('\n'.join(out))
  def check_slopes(self):
    for kq in range(4,14):
        for name in self.get_kq(kq):
          v=self[name]
          slopes=sign(diff(v))
          if sum(abs(diff(slopes)))>0:
             print name,' '.join(['%2d'%i for i in  slopes])
  def set_log(self):
    fig=pl.gcf()
    for ax in fig.axes:
        ax.set_xscale('log')
    pl.draw()
    return self
  def set_xlim(self,a,b):
    fig=pl.gcf()
    for ax in fig.axes:
        ax.set_xlim(a,b)
    pl.draw()
    return self
  def savefig(self):
    self.plot_betsqueeze()
    pl.savefig(self.filename.replace('.tfs','.png'))
    if len(self.get_acb(1))>0:
      self.plot_knobs()
      pl.savefig(self.filename.replace('.tfs','_knobs.png'))
    return self



