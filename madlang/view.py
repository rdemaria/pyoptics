"""
view
======

Provides:
  1. an high level interface for data object.
  2. prototype based object system
  3. support for property and autobind protocol

"""

def islist(i): return isinstance(i,list) or isinstance(i,tuple)
def isstr(i):  return hasattr(i,'startswith',)
def ispair(item):
  return islist(item) and len(item)==2 and isstr(item[0])
def isdict(i): return hasattr(i,'items')
def isview(i): return hasattr(i,'_data') and hasattr(i,'_proto')
def isobj(i):  return hasattr(i,'__dict__')


import sys
class EventHandler(object):
  """Dummy event handler"""
  def __init__(self):
    self.actors=[]
  def push(self,*args):
    for i in self.actors:
      try:
        i.run()
      except:
        print "Error in %s" % i
        print sys.exc_info()
  def register(self,actor):
    self.actors.append(actor)


class Actor(object):
  def __init__(self,message,method,target,active=True,args=()):
    self.active=active
    self.message=message
    self.target=target
    self.method=method
    self.arg=args
  def run(self,msg):
    if self.active and self.check(msg):
      getattr(self.target,self.method)(*args)
  def check(self,msg):
    res=True
    for i in range(len(self.message)):
      if msg[i]!=self.message[i]:
        res=False
        break
    return res
  def __repr__(self):
    return 'Actor(%s,%s,%s)' % (self.message,self.method,self.target)


class view(object):
  """ Generic data user interface"""
  _objectified=True
  _eventhandler=EventHandler()
  def __init__(self,data=None,proto=None):
    if data is None:
      data={}
    elif hasattr(data,'__dict__'):
      data=data.__dict__
    if proto is None:
      proto=[]
    self._data=data
    self._proto=proto
  def __getitem__(self,k):
    try:
      v=self._data[k]
    except KeyError:
      for p in self._proto:
        try:
          v=p[k]
          break
        except KeyError:
          pass
      else:
        v=self.__dict__[k]
    if hasattr(v,'_get'):
      return v._get(self,k)
    else:
      return v
  def __setitem__(self,k,v):
    if k in self._data and  hasattr(self._data[k],'_set'):
      self._data[k]._set(v,self,k)
    else:
      if hasattr(v,'_bind'):
        v=v._bind(self,k)
      self._data[k]=v
    self._eventhandler.push(self,'statechange')
  def __delitem__(self,k):
    if k in self._data:
      ov=self._data[k]
      if hasattr(ov,'_del'):
        ov._del(self,k)
      else:
        del self._data[k]
      self._eventhandler.push(self,'statechange')
    else:
      raise KeyError
  def __contains__(self,k):
    return k in self._data
  def items(self):
    if hasattr(self._data,'items'):
      return self._data.items()
    else:
      return []
  def _keys(self):
    return [ i[0] for i in self._data.items() if ispair(i) ]
  def _getAttributeNames(self):
    if self._objectified:
      out=self._keys()
      for i in self._proto:
        out.extend(i._getAttributeNames())
    else:
      out=[]
    out.extend(dir(self))
    return out
#  def __dir__(self):
#    return self._getAttributeNames()
  def __getattribute__(self,k):
    if k.startswith('_'):
      return object.__getattribute__(self,k)
    if object.__getattribute__(self,'_objectified'):
      try:
        return self[k]
      except KeyError:
        return object.__getattribute__(self,k)
    return object.__getattribute__(self,k)
  def __setattr__(self,k,v):
    if k.startswith('_'):
      object.__setattr__(self,k,v)
    elif self._objectified:
      self[k]=v
    else:
      object.__setattr__(self,k,v)
  def __delattr__(self,k):
    if k.startswith('_'):
      object.__delattr__(self,k,v)
    elif self._objectified:
      del self[k]
    else:
      object.__delattr__(self,k,v)
  def __repr__(self):
    """long representation"""
    attrs=[ k for k in self._keys() if not k.startswith('_')]
    if 'name' in attrs:
      name=' "%s"' % self['name']
    else:
      name=' %s' % id(self)
    if self._proto:
      proto='[%s]' % ','.join([str(i) for i in self._proto])
    else:
      proto=''
    out=['<%s%s%s' % (self.__class__.__name__, name,proto)]
    for k in attrs:
      v=self._data[k]
      if hasattr(v,'shape'):
        v='<array %s%s>' % (v.dtype.name,list(v.shape))
      out.append('  %-25s = %s' % (k,v))
    out[-1]=out[-1]+' >'
    if len(out)>120:
      out=out[:15]+['  ...']+out[-15:]
    return '\n'.join(out)
  def __str__(self):
    """short representation"""
    if 'name' in self._data:
      name=' "%s"' % self['name']
    else:
      name=' %s' % id(self)
    return '<%s%s>' % (self.__class__.__name__, name)

