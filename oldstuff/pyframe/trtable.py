import transport as _transport
from numpy import zeros,array,where,tan,ones
from matching import match
from utils import myflatten


_cols=_transport.colname().split()
_sindex=_cols.index('s')
_lindex=len(_cols)-_sindex

from opticstable import opticstable

class trtable(opticstable):

  def __init__(self,elems,full=True,inc_start=True,inc_stop=True):
    """ Class for optics computations. It allows to create a line of combined
    function magnets and evaluate the optics properties. The optics functions
    can be accessed by attribute name (e.g t.betx). Computation can be made
    (e.g. dqx=sum(t(kn1l*betx))/4/pi ). Element properties and function at
    the beginnig of the element can be accessed and modiefied (e.g t.mq.l=3).
    Selections can be used (e.g t[t.s>t.mq11.s].betx).

    elems
      is a string of space separated elements names. Element name can be
      any legal python keyowrds i.e [A-z][A-z0-9_]*

    full
      it create shortcut to access line properties by elems name. For line
      with a large number of names, the name access can create a performance
      slowdown.

    >>> from trtable import trtable
    >>> t=trtable('qf d qd d')
    >>> t.qf.kn1l=1
    >>> t.qd.kn1l=-1
    >>> t.d.l=1
    >>> t=t.track()
    """
    if not hasattr(elems,'__iter__'): elems=elems.split()
    if inc_start:
      elems.insert(0,'start')
    if inc_stop:
      elems.append('end')
    opticstable.__init__(self,elems,_cols,full=full)
    (self.br,self.gr,self.betx,self.bety)=(1,7.E6/938,1,1)

  def _get_array(self,start=None,stop=None):
    if start:
      start=where(self//start)[0][0]
    if stop:
      stop=where(self//start)[0][-1]+1
    return self._data[0][:,start:stop].T

  def track(self,stop=None,start=None):
    """ Track orbit, dispersion, beta function from initial conditions
    >>> from trtable import trtable
    >>> t=trtable('d')
    >>> t.d.l=23
    >>> t.start.betx=0.25
    >>> t.start.bety=0.55
    >>> t=t.track()
    >>> t.end.betx, t.end.bety
    (2116.25, 962.36818181818171)
    """
    # loop accept v(row,cols)
    _transport.track(self._get_array(start,stop))
    return self

  def ptrack(self,stop=None,start=None):
    _transport.ptrack(self._get_array(start,stop))
    return self

  def trace(self,stop=None,start=None):
    return _transport.trace(self._get_array(start,stop))

  def tmatrix(self,stop=None,start=None):
    """ Extract transfer matrix of the line. Optionally start and
      end element can be provided. Last element is not evaluated
    >>> from trtable import trtable
    >>> t=trtable('q e')
    >>> t.q.l=1
    >>> t.tmatrix()[0,1]==t.q.l
    True
    >>> from math import cos,cosh,sqrt
    >>> t.q.kn1l=1
    >>> t.tmatrix()[0,0]==cos(sqrt(t.q.kn1l/t.q.l)*t.q.l)
    True
    >>> t.tmatrix()[2,2]==cosh(sqrt(t.q.kn1l/t.q.l)*t.q.l)
    True
    """
    return _transport.tmatrix(self._get_array(start,stop))

  def split(self,n):
    oldnames=array(self._row_names[1:-1]) # remove start end
    r=len(oldnames)
    newnames=zeros(r*n,dtype=oldnames.dtype)
    for i in range(n):
      newnames[i::n]=oldnames
    t=trtable(newnames.tolist(),full=self._full)
    for c in myflatten(self._col_groups):
      for i in range(n):
        getattr(t,c)[i+1:-1:n]=getattr(self,c)[1:-1]
      getattr(t,c)[0]=getattr(self,c)[0]
      getattr(t,c)[-1]=getattr(self,c)[-1]
    t.l=t.l/n
    t.kn1l=t.kn1l/n
    t.ks1l=t.ks1l/n
    t.kn0l=t.kn0l/n
    t.ks0l=t.ks0l/n
    return t

  def __neg__(self):
    return reverse_table(self)

  def sbend_to_rbend(self,sel):
    """Add edgefocusing to selection
    >>> t=trtable('q b f f b q')
    >>> t.b.l=1; t.b.kn0l=1; t.l[2]=0
    >>> t2=t.sbend_to_rbend(t//'b')
    >>> t2.dump('.','l kn1l')
    name         l            kn1l        
    start           0.           0.       
    q               0.           0.       
    b_edge          0.           0.       
    b               0.           0.       
    b_edge          0.           0.       
    f               0.           0.       
    f               0.           0.       
    b_edge          0.         546.302e-03
    b               1.000        0.       
    b_edge          0.         546.302e-03
    q               0.           0.       
    end             0.           0.       
    """
    sel=(where(sel)[0]).tolist()
    names=self._row_names
    new_names=names[:]
    for i in sel:
      new_names[i]=[names[i]+'_edge',names[i],names[i]+'_edge']
    new_names=list(myflatten(new_names))
    edges=zeros(len(new_names),dtype=bool)
    sbends=zeros(len(new_names),dtype=bool)
    for i in range(len(sel)):
      edges [ sel[i]+i*2   ]= True
      sbends[ sel[i]+i*2+1 ]= True
      edges [ sel[i]+i*2+2 ]= True
    new_t=trtable(new_names,inc_start=False,inc_stop=False)
    for cname in myflatten(self._col_groups):
      col=getattr(new_t,cname)
      col[- edges]=getattr(self,cname)
    angle=new_t.kn0l[sbends]
    l=new_t.l[sbends]
    pos=l>0
    kick=( (angle/l)*tan(angle/2) )[pos]
    idx=where(edges)[0][0::2][pos]
    new_t.kn1l[idx]=kick
    new_t.kn1l[idx+2]=kick
    return new_t

  def reverse_trtable(t):
    names=t._row_names[:]
    names.pop()
    names.pop(0)
    names.reverse()
    tn=trtable(names)
  #  tn._all[:_sindex,:]=t._all[:_sindex,::-1]
    tn._all[:,:]=t._all[:,::-1]
    tn.alfx=-tn.alfx
    tn.alfy=-tn.alfy
    tn.dpx=-tn.dpx
    tn.dpy=-tn.dpy
    tn.px=-tn.px
    tn.py=-tn.py
    tn.s=-tn.s
    return tn

def from_tfstable(t,full=True,kicker=False):
  names=t._row_names
  l=trtable(names,full=full)
  if kicker:
    begcolm ='l k1l angle hkick vkick'.split()
    begcol  ='l kn1l kn0l kn0l ks0l'.split()
  else:
    begcolm ='l k1l angle'.split()
    begcol  ='l kn1l kn0l'.split()
  endcol='s betx bety alfx alfy dx dpx dy dpx x px y py mux muy'.split()
  for n,o in zip(begcol,begcolm):
    if hasattr(t,o):
      getattr(l,n)[1:-1]+=getattr(t,o)
      getattr(l,n)[0]+=getattr(t,o)[0]
      getattr(l,n)[-1]+=getattr(t,o)[-1]
  for n in endcol:
    if hasattr(t,n):
      getattr(l,n)[2:]=getattr(t,n)
      getattr(l,n)[0]=getattr(t,n)[0]
      getattr(l,n)[1]=getattr(t,n)[0]
#  l.add_col('name',array(names))
  return l




if __name__=='__main__':
  import doctest
  doctest.testmod()
  t=trtable('f b d d b f')
  t.l=1; t.start.l=0; t.end.l=0
  t.b.kn0l=1
  t.f.kn1l=1
  t.d.kn1l=-1
#  t.pl=t.track().plotbeta()
  from expr import expr
#  t.on_set(expr('track();pl.update();'))
  


