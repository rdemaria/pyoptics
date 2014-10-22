import cProfile,pstats
from tfstable import tfstable
from trtable import from_tfstable

def test():
  o=tfstable('twiss.lhcb1.data',rowattr=False)
  c=from_tfstable(o,rowattr=False)
  t=c.copy(rowattr=False)
#  return o,c,t

cProfile.run('test()','test.py.data')
p=pstats.Stats('test.py.data').strip_dirs()
p.sort_stats('cumulative').print_stats()
