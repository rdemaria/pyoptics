from tfsdata import tfsdata
from opticstable import opticstable
from plotfunc import plotfunc
from utils import pythonname as _pythonname
from trtable import trtable
from numpy import newaxis

class tfstable(opticstable,plotfunc):
  _is_s_begin=False
  def __init__(self,filename,full=True):
    tab=tfsdata(filename)
    col_names=tab.labels
    if 'name' in tab.labels:
      elems=map(_pythonname,tab.data['name'])
    else:
      elems=['l%05d' % i for i in range(tab.len) ]
    data={}
    opticstable.__init__(self,elems,full=full)
    for i in tab.data:
      self._add_cols([i],data=tab.data[i])
    if self._full:
       self._mkfull()

  def to_trtable(t,full=True,kicker=False):
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
  self=tfstable('twiss.ir5b1.data.gz')
  self.show('mqml','s bet')
  self.show((self.s > 4) & (self.s < 10),'s bet')
  print 'tfstable.py:  test OK'
