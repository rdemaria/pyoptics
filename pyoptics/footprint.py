from pydataobj import *
from tfsdata import *
from optics import errors_getmn,getmn
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.delaunay
import numpy as np
import os
from glob import glob
import re
import gzip

def find_res_xcross(m,n,q,xs,y1,y2,out):
#  print 'x=%d'% (xs),
  if n!=0:
    m,n,q,xs,y1,y2=map(float,(m,n,q,xs,y1,y2))
    ys=(q-m*xs)/n
#    print 'ys=%g'%ys,
#    print '%g<=ys<=%g'%(y1,y2),
    if ys>=y1 and ys<=y2:
#      print [xs,ys],
      out.append((xs,ys))
#  print

def find_res_ycross(m,n,q,ys,x1,x2,out):
#  print 'y=%d'% (ys),
  if m!=0:
    m,n,q,ys,y1,y2=map(float,(m,n,q,ys,x1,x2))
    xs=(q-n*ys)/m
#    print 'xs=%g'%xs,
#    print '%g<=xs<=%g'%(x1,x2),
    if xs>=x1 and xs<=x2:
#      print [xs,ys]
      out.append((xs,ys))
#  print

def get_res_box(m,n,l=0,qz=0,a=0,b=1,c=0,d=1):
  """get (x,y) coordinates of resonance lines with
     m,n,q:   resonance integers with m*qx+n*qy=q
     l,qz:    order l of resonance sideband with frequency qz
     a,b,c,d: box parameters=tune range, 
              explicitly a<qx<b and c<qy<d 
  """
  order=int(np.ceil(abs(m)*max(abs(a),abs(b))+abs(n)*max(abs(c),abs(d))))
  out=[]
  for q in range(-order,+order+1):
    q=q-l*qz
    points=[]
#    print '%2d*Qx+%2d*Qy=%2d' % (m,n,q)
    find_res_xcross(m,n,q,a,c,d,points)#find endpoint of line (a,ys) with c<ys<d
    find_res_xcross(m,n,q,b,c,d,points)#find endpoint of line (b,ys) with c<ys<d
    find_res_ycross(m,n,q,c,a,b,points)#find endpoint of line (xs,c) with a<xs<b
    find_res_ycross(m,n,q,d,a,b,points)#find endpoint of line (xs,d) with a<xs<b
    points=list(set(points))
    if len(points)>1:
      out.append(points)
  return out

def plot_res_box(m,n,l=0,qz=0,a=0,b=1,c=0,d=1,color='b',linestyle='-'):
  points=get_res_box(m,n,l,qz,a,b,c,d)
  for c in points:
    x,y=zip(*c)
    pl.plot(x,y,color=color,linestyle=linestyle)


def plot_res_order_box(o,l=0,qz=0,a=0,b=1,c=0,d=1,c1='b',lst1='-',c2='b',lst2='--',c3='g'):
  """plot resonance lines for order o and sidbands of order l and frequency qz"""
  for m,n in getmn(o,'b'):
    # print 'b%s: m=%d n=%d'%(o,m,n)
    plot_res_box(m,n,l=0,qz=0,a=a,b=b,c=c,d=d,color=c1,linestyle=lst1)
    if(l!=0):#sidebands
      for ll in +abs(l),-abs(l):
        plot_res_box(m,n,l=ll,qz=qz,a=a,b=b,c=c,d=d,color=c3,linestyle=lst1) 
  for m,n in getmn(o,'a'):
    # print 'a%s: m=%d n=%d'%(o,m,n)
    plot_res_box(m,n,l=0,qz=0,a=a,b=b,c=c,d=d,color=c2,linestyle=lst2)
    if(l!=0):#sidebands
      for ll in +abs(l),-abs(l):
        plot_res_box(m,n,l=ll,qz=qz,a=a,b=b,c=c,d=d,color=c3,linestyle=lst2) 


def plot_res(m,n,l=0,qz=0,color='b',linestyle='-'):
  a,b=pl.xlim()
  c,d=pl.ylim()
  points=get_res_box(m,n,l,qz,a,b,c,d)
  for c in points:
    x,y=zip(*c)
    pl.plot(x,y,color=color,linestyle=linestyle)
  pl.xlim(a,b)
  pl.ylim(c,d)
def plot_res_order(o,l=0,qz=0,c1='b',lst1='-',c2='b',lst2='--',c3='g'):
  a,b=pl.xlim()
  c,d=pl.ylim()
  plot_res_order_box(o,l,qz,a,b,c,d,c1,lst1,c2,lst2,c3)
  pl.xlim(a,b)
  pl.ylim(c,d)

#a=c=0.28;b=d=0.34;
#map(plot_res_order,range(1,15))




def mkranges(nsigmax=12,nangles=7):
  ranges=[]
  for i in range(nangles):
    ranges.append([0,i])
  for i in range(nangles):
    ranges.append(slice(1+i,nangles*nsigmax+1,nangles))
  for i in range(nsigmax):
    ranges.append(slice(1+nangles*i,1+nangles*(i+1)))
  return ranges

def plot_grid(t,nsigmax=12,nangles=7,lw=1):
  ranges=mkranges(nsigmax,nangles)
  for i in ranges:
    if hasattr(i,'step'):
      lw= i.step==nangles and i.start/2. or 1
    pl.plot(t.x[i],t.y[i],'-k',lw=lw)

def plot_footprint(t,name='',nsigmax=6,nangles=7,wp=None,spread=0.01,color=None):
  ranges=mkranges(nsigmax,nangles)
  lw=1
  out=[]
  lbl=True
  if color==None: color=colorrotate()
  for i in ranges:
    if lbl:
       p=pl.plot(t.tunx[i],t.tuny[i],'-%s'%color,lw=lw,label=name)
       lbl=False
    else:
       p=pl.plot(t.tunx[i],t.tuny[i],'-%s'%color,lw=lw)
    out.append(p[0])
  pl.ylabel('$Q_y$')
  pl.xlabel('$Q_x$')
  pl.grid(True)
  if wp!=None:
    qx,qy=wp
    pl.xlim(qx-spread,qx+spread)
    pl.ylim(qy-spread,qy+spread)
  return out

mycolors=list('rgbcm')

def colorrotate():
  c=mycolors.pop(0);mycolors.append(c)
  return c


def mkplot(l,n,s):
  l,n,s=map(int,(l,n,s))
  fn='jobslhc3_%d_mb_%d_v2/dynaptune_%d.tfs'% (l,n,s)
  t=dataobj(tfsdata.open(fn))
  c=mycolors.pop(0);mycolors.append(c)
  plot_footprint(t,fn,nsigmax=12,color='c')
  return t

def triangulate(t):
  tr=matplotlib.delaunay.triangulate.Triangulation(t.tunx,t.tuny)
  for i in tr.triangle_nodes:
    plot(t.tunx[i],t.tuny[i])

#*********************** START: functions to plot resonances in xy space - takes as input dynaptune files *************************
def get_res(o,l,qz,tunexmin,tunexmax,tuneymin,tuneymax):
  '''return the resonances lines (skew and normal -> kind='ab') which cross the area [(tunexmin,tunexmax),(tuneymin,tuneymax)]'''
  res_list=[]
  for m,n in getmn(o,'ab'):#get m,n for resonance of order o
    rbox=get_res_box(m,n,l,qz,tunexmin,tunexmax,tuneymin,tuneymax)
    if(len(rbox)>0):
      for qq in rbox:
        q=m*qq[0][0]+n*qq[0][1]+l*qz
        res_list.append((o,m,n,int(round(q)),l,qz))
  return res_list

def get_xy_dat_res(data,res,eps):
  '''get the particles with a abs(m*tunex+n*tuney-r)=delta<eps'''
  (o,m,n,q,l,qz)=res
  idx=(abs(m*(data.tunx)+n*(data.tuny)+l*qz-q)<eps)
  xres=data.x[idx]
  yres=data.y[idx]
  tunxres=data.tunx[idx]
  tunyres=data.tuny[idx]
  ftype=[('o',int),('m',int),('n',int),('q',int),('l',int),('qz',float),('eps',float),('x',float),('y',float),('tunx',float),('tuny',float),('delta',float)]
  xyres=np.zeros(len(xres),dtype=ftype)
  xyres['o']=o
  xyres['m']=m
  xyres['n']=n
  xyres['l']=l
  xyres['q']=q
  xyres['qz']=qz
  xyres['eps']=eps
  xyres['x']=xres
  xyres['y']=yres
  xyres['tunx']=tunxres
  xyres['tuny']=tunyres
  xyres['delta']=abs(m*tunxres+n*tunyres+l*qz-q)
  return xyres

def plot_res_xy(datan,omin=1,omax=10,l=0,qz=0,eps=1.e-2,plxlim=None,plylim=None):
  '''plot resonances in x,y space up to order omax and with maximum 
  distance abs(m*tunex+n*tuney+l*qz-r)=delta<eps. The data should be given
  in number of sigma -> x[sigma], y[sigma]'''
  fig=pl.gcf()
  gs=gridspec.GridSpec(1,1)
  gs.update(hspace=0.4,wspace=0.5)
  ax=fig.add_subplot(gs[0])
  cmin=0
  cmax=1
  cb=None
  print omax
  for reso in range(omin,int(round(omax+1))):
    lres=get_res(reso,l,qz,min(datan.tunx),max(datan.tunx),min(datan.tuny),max(datan.tuny))
    for ll in lres:
      xyres=get_xy_dat_res(datan,ll,eps)
      if(len(xyres['x'])>0):
        print ll
        if cb is None:
          s = ax.scatter(xyres['x'],xyres['y'],c=(1-xyres['delta']/eps),marker='.',edgecolors='none')
          s.set_clim([cmin,cmax])
          cb=fig.colorbar(s)
        else:
          s = ax.scatter(xyres['x'],xyres['y'],c=(1-xyres['delta']/eps),marker='.',edgecolors='none')
          s.set_clim([cmin,cmax])
  #      pl.plot(xyres['x'],xyres['y'],marker='.',linestyle='None',color='b')
  ax.tick_params(axis='x', pad=5)
  ax.tick_params(axis='y', pad=5)
  pl.xlabel(r'$\sigma_x$',fontsize=16,labelpad=14)
  pl.ylabel(r'$\sigma_y$',fontsize=16,labelpad=14)
  pl.xlim(plxlim)
  pl.ylim(plylim)
  cb.set_label(r'$1-\Delta$res/eps with eps=%s'%round(eps,4))

#*********************** END: functions to plot resonances in xy space - takes as input dynaptune files *************************

#fn='/afs/cern.ch/project/uslarp/rdemaria/work/slhc3_1/dynamptune_thin_8282.tfs'
#t=dataobj(tfsdata.open(fn))
#c=mycolors.pop(0);mycolors.append(c)
##plot_grid(t,12)
#plot_footprint(t,fn,nsigmax=6,c='k')
#
#
#fn='/afs/cern.ch/project/uslarp/rdemaria/work/slhc3_2/dynamptune_thin_8282.tfs'
#t=dataobj(tfsdata.open(fn))
##plot_grid(t,12)
#plot_footprint(t,fn,nsigmax=6,c='k')
#
#
#fn='/afs/cern.ch/project/uslarp/rdemaria/work/slhc3_3/dynamptune_thin_8282.tfs'
#t=dataobj(tfsdata.open(fn))
#c=mycolors.pop(0);mycolors.append(c)
##plot_grid(t,12)
#plot_footprint(t,fn,nsigmax=6,c='k')
#
#mkplot(1,8282,1);mkplot(2,8282,1); mkplot(3,8282,1);


#for i in range(1,61):
#  try:
#    mkplot(3,8282,i)
#  except:
#    print 'error on',i
#
#xlim(0.298,0.312)
#ylim(0.3165,0.3205)

def change_alpha(ll,alpha):
  for a in ll:
    hasasttr(a,'set_alpha') and a.set_alpha(alpha) or change_alpha(a,alpha)

#regexp=re.compile('dynaptune_(..)(..)_(....)_dqx(......)_dqy(......)(.?)')
#def mkfn(bb,ir,squeeze,dqx,dqy,opt):
#  return 'dynaptune_%s%s_%s_dqx%+6.3f_dqy%+6.3f%s'%(bb,ir,squeeze,dqx,dqy,opt)
#
#def load(bb,ir,squeeze,dqx,dqy,opt):
#  t=dataobj(tfsdata.open(mkfn(bb,ir,squeeze,dqx,dqy,opt)))
#  t.tunx-=dqx; t.tuny-=dqy
#  return t
#
#data={}
#for i in glob('dynaptune*'):
#  bb,ir,squeeze,dqx,dqy,opt=regexp.match(i).groups()
#  dqx=float(dqx)
#  dqy=float(dqy)
#  print bb,ir,squeeze,dqx,dqy,opt
#  t=load(bb,ir,squeeze,dqx,dqy,opt)
#  data[(bb,ir,squeeze,opt)]=t


