import matplotlib.pyplot as pl


from optics import optics

class IRPlots(object):
  def __init__(self,name):
    self.name=name
    self.plot()
  def plot(self):
    self.p1=self.plotn(1)
    self.p2=self.plotn(2)
    return self
  def plot(self):
    self.p1=self.plotn(1)
    self.p2=self.plotn(2)
    return self
  def plotn(self,n):
    name="%sb%d"%(self.name,n)
    pl.figure(name)
    p=optics.open("twiss_%s.tfs"%name).plotbeta(newfig=False)
    pl.title(name)
    p.wx_autoupdate()
    return p
  def save(self,figname):
    f1='%sb1_%s.png'%(self.name,figname)
    f2='%sb2_%s.png'%(self.name,figname)
    self.p1.figure.savefig(f1)
    self.p2.figure.savefig(f2)
    print f1,f2
  def scalebeta(self,betamax):
    self.p1.left.set_ylim((0,betamax))
    self.p2.left.set_ylim((0,betamax))
    self.p1.run()
    self.p2.run()




