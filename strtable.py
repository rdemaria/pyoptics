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
    return self.get_vars(r'mu[xy]ip[1-8]b[12]_[lr]')
  def plot_triplet(self,n1,n2,x=None):
    scale=StrTable.scale
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    ks=self.get_vars('kqx.*r')
    for k in ks:
      pl.plot(xv,abs(self[k][n1:n2]*scale),label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel('k')
  def plot_2in1(self,kq,n1,n2,x=None,sign=False,ylab='k [T/m]'):
    scale=StrTable.scale
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
    for k in self.get_phases():
      kv=self[k][n1:n2]
      pl.plot(xv,kv,label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel('mu')
  def plot_squeeze(self,n1=0,n2=None,x=None):
    fig=pl.figure('squeeze',figsize=(16,12))
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
    self.plot_ipbeta(n1,n2,x=x)
    return self
  def plot_betsqueeze(self,val=0.44):
    x=self.get_vars('betxip')[0]
    n1=where(self[x]==val)[0][-1]
    n2=None
    self.plot_squeeze(n1,n2,x)
    return self
  def linint(t,name,val1,val2,x,par='ttt'):
    n1=where(t[x]==val1)[0][-1]
    n2=where(t[x]==val2)[0][-1]
    tmp= "%s:=(%g)*(1-%s)+(%g)*(%s);"
    print tmp%(name,t[name][n1],par,t[name][n2],par)
  def button_press(self,event):
    self.event=event
  def poly_fit(self,var,order,param="betxip8b1",slope0=[]):
    scale=StrTable.scale
    x=self[param]; y=self[var]
    x0=[x[0],x[-1]]
    y0=[y[0],y[-1]]
    xp0=[]; yp0=[]
    for idx in slope0:
        xp0.append(x[idx])
        yp0.append(0)
    pol=poly_fit(order,x,y,x0,y0,xp0,yp0)
    out="%s:=%s;"%(var, poly_print(pol,x=param,power='^'))
    n=re.match('[kqtxl]+([0-9]+)\.',var)
    if n is None:
        n=3
    else:
        n=int(n.groups()[0])
    pl.subplot(3,4,n-2)
    yv=poly_val(pol,x)
    print ' '.join(['%2d'%i for i in  sign(diff(yv))])
    if n<11:
        yv=abs(yv)
    pl.plot(yv*scale)
    print out
  def check_slopes(self):
    for kq in range(4,14):
        for name in self.get_kq(kq):
          v=self[name]
          slopes=sign(diff(v))
          if sum(abs(diff(slopes)))>0:
             print name,' '.join(['%2d'%i for i in  slopes])




