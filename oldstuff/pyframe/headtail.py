from view import view
from numpy import fromfile


c=slice(None)
class hdtltable(view):
  _is_s_begin=False
  def __init__(self,filename,slices,turns=None):
    t=fromfile(filename,sep='\n')
    if not turns:
      turns=t.size/6/slices
    t=t.reshape( (turns,slices,6))
    view.__init__(self,
        all=t,
        z=(t,(c,c,0)),
        x=(t,(c,c,1)),
        y=(t,(c,c,2)),
        sx=(t,(c,c,3)),
        sy=(t,(c,c,4)),
        d=(t,(c,c,5))
        )

if __name__=='__main__':
  print 'tfstable.py:  test OK'
