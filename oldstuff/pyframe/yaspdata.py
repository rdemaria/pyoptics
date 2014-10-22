def myopen(fn):
  if fn.endswith('.gz'):
    import gzip
    return gzip.open(fn)
  else:
    return open(fn)


import numpy as n

yasptypes={'%s': str, '%d': int, '%f': float}
coltypes={'BEAM': 'i',
 'HW-STATUS': 'i',
 'KICK': 'd',
 'NAME': 'S20',
 'PLANE': 'S1',
 'POS': 'd',
 'RMS': 'd',
 'STATUS': 'i',
 'STATUS-TAG': 'S20',
 'STRENGTH-NAME': 'S20',
 'SUM': 'd'}

def readdata(o,cols,count):
  ds=[]
  for i in range(count):
    ds.append(o.readline().split())
  ds=zip(*ds)
  d={}
  for name,data in zip(cols,ds):
    d[name.lower()]=n.array(data,dtype=coltypes[name])
  return d

class yaspdata(object):
  def __init__(self,filename):
    self.filename=filename
    o=myopen(filename)
    c=''
    param={}
    dataset={}
    l=o.readline()
    while l:
      c=l[0]
      if c=='@':
        name,name,type,data=l.strip().split(None,4)
        param[name]=yasptypes[type](data)
      elif c=='#':
        name,name,type,data=l.strip().split(None,4)
        cols=o.readline().strip().split();cols.pop(0)
        if name=='MONITOR':
          dataset['monitor-h']=readdata(o,cols,param['MONITOR-H-NUM'])
          dataset['monitor-v']=readdata(o,cols,param['MONITOR-V-NUM'])
        elif name=='CORRECTOR':
          dataset['corrector-h']=readdata(o,cols,param['CORRECTOR-H-NUM'])
          dataset['corrector-v']=readdata(o,cols,param['CORRECTOR-V-NUM'])
      l=o.readline()
    self.dataset=dataset
    self.param=param








