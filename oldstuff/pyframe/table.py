from frame import frame
from view import view
from utils import mystr as _mystr
from utils import myflatten
from numpy import zeros,array, ones, where, r_
import re



def _mkidx(lst):
  """ Create a dictionary of elements and element index in the list
  """
  d={}
  for i,n in enumerate(lst):
    d.setdefault(n,[]).append(i)
  return d

class stable(frame):
  _name_char=12
  _entry_char=12
  _entry_prec=3
  """ Container class where attributes are organized in row and colums.
      Interfaces are:
      - _row_names : list of row names
      - _col_names : list of column names
      Provides:
      - pattern : create an index of rows mathing the pattern
  """
  def __init__(self,row_names=[],col_names=[],full=True):
    """
    >>> t=stable()
    >>> t.show()
    name        
    >>> t.add_cols('c1 c2 c3')
    >>> t.show()
    name         c1           c2           c3          
    >>> t.add_rows('r1 r2 r3')
    >>> t.show()
    name         c1           c2           c3          
    r1              0.           0.           0.       
    r2              0.           0.           0.       
    r3              0.           0.           0.       
    >>> t.c1=1; t.r2=2; t.r3.c2=3
    >>> t.show()
    name         c1           c2           c3          
    r1              1.000        0.           0.       
    r2              2.000        2.000        2.000    
    r3              1.000        3.000        0.       
    """
    frame.__init__(self)
    self._row_names=[]
    self._row_index={}
    self._data=[]
    self._col_groups=[]
    self._full=full
    if row_names:
      self.add_rows(row_names)
    if col_names:
      self.add_cols(col_names)

  def add_cols(self,col_names,data=None,dtype=float):
    self._add_cols(col_names,data=data,dtype=dtype)
    if self._full:
      self._mkfull()

  def _add_cols(self,col_names,data=None,dtype=float):
    if not hasattr(col_names,'__iter__'):  col_names=col_names.split()
    if data is None:
      ar=zeros( ( len(self._row_names), len(col_names) ),dtype=dtype )
    else:
      ar=data
    self._data.append(ar)
    for j,n in enumerate(col_names):
      if len(ar.shape)==1:
        idx=slice(None)
      else:
        idx=(slice(None),j)
      self._overattr(n,view(ar, idx))
    self._col_groups.append(col_names)

  def add_rows(self,row_names,idx=None):
    if not hasattr(row_names,'__iter__'):  row_names=row_names.split()
    if idx is None:
      idx=len(self)
    lrow=len(row_names)
    if lrow>0:
      self._row_names.insert(idx,row_names)
      self._row_names=list(myflatten(self._row_names))
    new_data=[]
    for d in self._data:
      newd=zeros( (lrow,d.shape[1]), dtype=d.dtype)
      new_data.append(r_[ d[:idx], newd,  d[idx:] ])
    self._data=new_data
    for cnames,d in zip(self._col_groups,self._data):
      for cname in cnames:
        self.__dict__[cname].data=d
    self._row_index=_mkidx(self._row_names)
    if self._full:
      self._mkfull()

  def _mkfull(self):
    bclass=self.__class__.__base__
    for rname,idx in self._row_index.items():
      if len(idx)==1:
        idx=idx[0]
      rclass=bclass()
      for cnames,d in zip(self._col_groups,self._data):
        for cidx,cname in enumerate(cnames):
          if len(d.shape)==1:
            rclass._overattr(cname,view(d,idx) )
          else:
            rclass._overattr(cname,view(d,(idx,cidx)) )
      self._overattr(rname,rclass)

  def __len__(self):
    return len(self._row_names)

  def _extract(self,k):
    t=table(full=self._full)

  def __getitem__(self,k):
    if isinstance(k,slice):
      return k
    else:
      return frame.__getitem__(self,k)

#  def __repr__(self):
#    return self.dump()

#  __str__=__repr__

  def pattern(self,regexp):
    """ extract mask from pattern
    """
    c=re.compile(regexp)
    out=[c.search(n) is not None for i,n in enumerate(self._row_names)]
    return array(out)
  __floordiv__=pattern

  def dump(self,rows=None,cols=None):
    """show table
    """
    if rows is None:
      rows=ones(len(self._row_names),dtype=bool)
    elif isinstance(rows,str):
      rows=self.pattern(rows)
    rows=where(rows)[0]

    if cols is None:
      colsn=list(myflatten(self._col_groups))
      cols=[getattr(self,n) for n in colsn]
    if isinstance(cols,str):
      colsn=cols.split()
      cols=[self._eval(n) for n in cols.split()]

    out=[]
    rowfmt=['%%-%d.%ds' % (self._name_char,self._name_char)]
    rowfmt+=['%%-%d.%ds' % (self._entry_char,self._name_char)] * len(colsn)
    rowfmt=' '.join(rowfmt)
    out.append(rowfmt % tuple(['name'] + colsn  ) )
    if len(rows)==1:
      v=[ self._row_names[0] ]+ [ _mystr(c,self._entry_char) for c in cols ]
      out.append(rowfmt %  tuple(v))
    else:
      for i in rows:
        v=[ self._row_names[i] ]+ [ _mystr(c[i],self._entry_char) for c in cols ]
        out.append(rowfmt %  tuple(v))
    return '\n'.join(out)

  def show(self,rows=None,cols=None):
    print self.dump(rows=rows,cols=cols)

  def dumpraw(self,rows=None,cols=None):
    out=['::\n']
    for i in self.dump(rows=rows,cols=cols).split('\n'):
      out.append('  %s' % i)
    return '\n'.join(out)

if __name__=='__main__':
    import doctest
    doctest.testmod()
    t=stable()
    t.add_cols('c1 c2 c3')
    t.add_rows('r1 r2 r3')
    t.c1=1; t.r2=2; t.r3.c2=3



