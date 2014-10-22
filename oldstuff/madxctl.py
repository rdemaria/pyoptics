import os, pexpect
import sys
import re

class madxctl(object):
  def __init__(self,debug=True):
    self.child=pexpect.spawn('madx',maxread=200000,searchwindowsize=30)
    self.child.setecho(False)
    self.child.timeout=60
    if debug:
      self.child.logfile=sys.stdout
    self.prompt=['X: ==>\r\n', pexpect.TIMEOUT,'yes>:\r\n\r\n']
    self.send('')
    self.send('option,-echo;')
  def send(self,cmds):
    cmds=cmds.split(';')
    for cmd in cmds:
      self.child.sendline(cmd+';')
      self.child.expect(self.prompt)
      if self.child.after=='yes>:\r\n\r\n':
        self.noplot()
        self.prompt=['X: ==>\r\n', pexpect.TIMEOUT]
    return self.child.before
  def call(self,s):
    return self.send('call, file="%s";'%s)
  def noplot(self):
    self.send('no\nno\n')
  def quit(self):
    self.child.sendeof()
  def value(self,expr):
    b=self.send('value, %s;' % expr)
    p=madxParseValue(b)
    p=map(dict,p)
    p=[ (i['name'],eval(i['value'])) for i in p]
    return dict(p)
  def show(self,name):
    self.send('show, %s;' % name)
    return self.child.before
  def twiss_opt(self):
    self.send('select,flag=twiss,clear;')
    self.send('select, flag=twiss, column=name,parent,s,l,angle,k1l,betx,bety,alfx,alfy,dx,dpx,mux,muy,x,y,px,py,dy,dpy;')
  def twissdata(self,location,data='betx bety'):
    out=dict(location=location)
    for name in data.split():
      out[name]=self.value("table(twiss,%s,%s)" % (location,name))
    return out
  def twisstable(self):
    self.send('write,table=twiss;')
    return self.child.before

def madxParseValue(s):
  try:
    s=re.sub('\s','',s)
    out=[]
    regex=re.compile('([^:=]*)(:?=)(.*)')
    var='name eq value'
    for i in s.split(';'):
      if i:
        out.append(zip(var.split(),re.match(regex,i).groups()))
    return out
  except:
    raise 'Parse Error:','%s' % s


