import pylab as _p
import numpy as _n

lglabel={
    'betx':    r'$\beta_x$',
    'bety':    r'$\beta_y$',
    'dx':    r'$D_x [m]$',
    'dy':    r'$D_y [m]$',
    }

axlabel={
    's':       r'$s [m]$',
    'betx':    r'$\beta [m]$',
    'bety':    r'$\beta [m]$',
    'dx':    r'$D [m]$',
    'dy':    r'$D [m]$',
    'x':    r'$co [m]$',
    'y':    r'$co [m]$',
    }


def _mylbl(d,x): return d.get(x,r'$%s$'%x)

class qdplot(object):
  def __init__(self,t,x,yl,yr,idx,lattice,newfig,clist,pre):
    self.color={}
    self.left=None
    self.right=None
    self.lattice=None
    self.pre=None
    self.t,self.x,self.yl,self.yr,self.idx,self.clist=t,x,yl,yr,idx,clist
    self.xaxis=getattr(self.t,self.x)
    for i in self.yl+self.yr:
      self.color[i]=self.clist.pop(0)
      self.clist.append(self.color[i])
    if newfig:
      self.figure=_p.figure()
    else:
      self.figure=_p.gcf()
    self.figure.subplots_adjust(right=0.72)
    if lattice:
      self.lattice=self._new_axes()
#      self.lattice.set_autoscale_on(False)
      self.lattice.yaxis.set_visible(False)
    if yl:
      self.left=self._new_axes()
#      self.left.set_autoscale_on(False)
      self.left.yaxis.set_label_position('left')
      self.left.yaxis.set_ticks_position('left')
    if yr:
      self.right=self._new_axes()
#      self.right.set_autoscale_on(False)
      self.right.yaxis.set_label_position('right')
      self.right.yaxis.set_ticks_position('right')
    self.update()
  def _new_axes(self):
    if self.figure.axes:
      ax=self.figure.axes[-1]
      out=self.figure.add_axes(ax.get_position(),
          sharex=ax, frameon=False)
    else :
      out=self.figure.add_axes([.1,.1,.65,.8])
    return out
  def __repr__(self):
    return object.__repr__(self)
  def update(self):
    if callable(self.pre):
      self.pre()
    _p.ioff()
    self.lines=[]
    self.legends=[]
#    self.figure.lines=[]
#    self.figure.patches=[]
#    self.figure.texts=[]
#    self.figure.images = []
    self.figure.legends = []

    if self.lattice:
      self.lattice.patches=[]
      self._lattice(['k0l','kn0l','angle'],"#a0ffa0",'Bend h')
      self._lattice(['ks0l'],"#ffa0a0",'Bend v')
      self._lattice(['kn1l','k1l'],"#a0a0ff",'Quad')
      self._lattice(['hkick'],"#e0a0e0",'Kick h')
      self._lattice(['vkick'],"#a0e0e0",'Kick v')
    if self.left:
      self.left.lines=[]
      for i in self.yl:
        self._column(i,self.left,self.color[i])
    if self.right:
      self.right.lines=[]
      for i in self.yr:
        self._column(i,self.right,self.color[i])
    _p.xlabel(_mylbl(axlabel,self.x))
    self.figure.gca().set_xlim(min(self.xaxis[self.idx]),max(self.xaxis[self.idx]))
    _p.figlegend(self.lines,self.legends,'upper right')
    _p.grid(True)
#    self.figure.canvas.mpl_connect('button_release_event',self.button_press)
#    self.figure.canvas.mpl_connect('pick_event',self.pick)
    _p.ion()
    _p.draw()


#  def pick(self,event):
#    self.pickpos=_n.array([event.mouseevent.x,event.mouseevent.y])
#    self.pickname=event.artist.elemname
#    self.pickprop=event.artist.elemprop
#    print '\n',event.artist.elemname, event.artist.elemprop,

#  def button_press(self,mouseevent):
#    rel=_n.array([mouseevent.x,mouseevent.y])
#    dx,dy=self.pickpos/rel
#    print 'release'
#    self.t[self.pickname][self.pickprop]*=dy
#    self.t.track()
#    self.update()

  def _lattice(self,names,color,lbl):
    vd=0
    sp=self.lattice
    s=self.xaxis
    for i in names:
      if hasattr(self.t,i):
        vdname=i
        vd=getattr(self.t,i)[self.idx]+vd
    if vd is not 0:
      m=_n.abs(vd).max()
      if m>1E-10:
        c=_n.where(abs(vd) > m*1E-4)[0]
        if len(c)>0:
          if _n.all(self.t.l[c]>0):
            vd[c]=vd[c]/self.t.l[c]
            m=abs(vd[c]).max()
          vd[c]/=m
          if self.t._is_s_begin:
            plt=self.lattice.bar(s[c],vd[c],self.t.l[c],color='#aaffaa',picker=True)
          else:
            plt=self.lattice.bar(s[c]-self.t.l[c],vd[c],self.t.l[c],color='#aaffaa',picker=True)
          _p.setp(plt,facecolor=color,edgecolor=color)
          if plt:
            self.lines.append(plt[0])
            self.legends.append(lbl)
          for r,i in zip(plt,c):
            r.elemname=self.t._row_names[i]
            r.elemprop=vdname
        self.lattice.set_ylim(-1.5,1.5)

  def _column(self,name,sp,color):
    fig,s=self.figure,self.xaxis
    y=getattr(self.t,name)[self.idx]
    bxp,=sp.plot(s,y,color,label=_mylbl(lglabel,name))
    sp.set_ylabel(_mylbl(axlabel,name))
    self.lines.append(bxp)
    self.legends.append(_mylbl(lglabel,name))
    sp.autoscale_view()


class plotfunc(object):
  _is_s_begin=True

  def plot(self,yl='',yr='',x='s',idx=slice(None),clist='k r b g c m',lattice=True,newfig=True,pre=None):
    yl,yr,clist=map(str.split,(yl,yr,clist))
    out=qdplot(self,x,yl,yr,idx,lattice,newfig,clist,pre)
    return out
  def plotbeta(self,newfig=True):
    return self.plot('betx bety','dx dy',newfig=newfig)

  def plotcross(self,newfig=True):
    return self.plot('x y','dx dy',newfig=newfig)

  def plottune(self,lbl=''):
    _p.title(r"${\rm Tune} \quad {\rm vs} \delta$")
    _p.xlabel("$\delta$")
    _p.ylabel("Fractional tune")
    tt=r'$%s \rm{%s}$'
    _p.plot(self.deltap,self.q1-self.q1.round(),label=tt %('Q_x',lbl))
    _p.plot(self.deltap,self.q2-self.q2.round(),label=tt %('Q_y',lbl))
    qx=(self.q1-self.q1.round())[abs(self.deltap)<1E-15][0]
    qy=(self.q2-self.q2.round())[abs(self.deltap)<1E-15][0]
    _p.text(0.0,qx,r"$Q_x$")
    _p.text(0.0,qy,r"$Q_y$")
    _p.grid(True)
    _p.legend()

  def plotbetabeat(self,t1,dp='0.0003'):
    _p.title(r"$\rm{Beta beat: 1 - \beta(\delta=%s)/\beta(\delta=0)}$" % dp)
    _p.ylabel(r"$\Delta\beta/\beta$")
    _p.xlabel(r"$s$")
    _p.plot(self.s,1-t1.betx/self.betx,label=r'$\Delta\beta_x/\beta_x$')
    _p.plot(self.s,1-t1.bety/self.bety,label=r'$\Delta\beta_y/\beta_y$')
    _p.grid(True)
    _p.legend()

  def plotw(self,lbl=''):
    title(r"Chromatic function: %s"%lbl)
  #  ylabel(r"$w=(\Delta\beta/\beta)/\delta$")
    ylabel(r"$w$")
    xlabel(r"$s$")
    plot(self.s,self.wx,label=r'$w_x$')
    plot(self.s,self.wy,label=r'$w_y$')
    grid(True)
    legend()

