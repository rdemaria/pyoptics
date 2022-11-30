import os

import matplotlib.pyplot as pl
import matplotlib


from .optics import optics


class IRPlots(object):
    def __init__(self, name, basedir="."):
        self.name = name
        self.basedir = basedir
        self.plot()

    def plot(self):
        self.t1 = self.plotn(1)
        self.t2 = self.plotn(2)
        return self

    def plotn(self, n, lbl=""):
        name = "%sb%d" % (self.name, n)
        pl.figure(name)
        fn = os.path.join(self.basedir, "twiss_%s.tfs" % name)
        t = optics.open(fn).plotbeta(newfig=False)
        pl.title(name + lbl)
        if "wx" in matplotlib.get_backend().lower():
            t._plot.wx_autoupdate()
        else:
            t._plot.ani_autoupdate()
        return t

    def save(self, figname):
        f1 = "%sb1_%s.png" % (self.name, figname)
        f2 = "%sb2_%s.png" % (self.name, figname)
        self.t1._plot.figure.savefig(f1)
        self.t2._plot.figure.savefig(f2)
        print(f1, f2)
        return self

    def scalebeta(self, betamax):
        self.t1._plot.left.set_ylim((0, betamax))
        self.t2._plot.left.set_ylim((0, betamax))
        self.t1._plot.run()
        self.t2._plot.run()
        return self

    def scaledisp(self, dmin, dmax):
        self.t1._plot.right.set_ylim((dmin, dmax))
        self.t2._plot.right.set_ylim((dmin, dmax))
        self.t1._plot.run()
        self.t2._plot.run()
        return self

    def plotap(self, nlim=30, ref=7):
        irn = self.name[-1]
        for beam in "12":
            apfn1 = os.path.join(self.basedir, "ap_%sb%s.tfs" % (self.name, beam))
            fn = os.path.join(self.basedir, "twiss_%sb%s.tfs" % (self.name, beam))
            ap1 = optics.open(apfn1)
            t = optics.open(fn).plotap(ap1, nlim=nlim, ref=ref)
            pl.ylim(6.3, nlim)
            pl.title("ir%sb%s" % (irn, beam))
            pl.savefig(apfn1.replace(".tfs", ".png"))
        return self

    def plotcross(self, lss=True):
        irn = self.name[-1]
        for beam in "12":
            if lss:
                n1 = "e.ds.l%s.b%s" % (irn, beam)
                n2 = "s.ds.r%s.b%s" % (irn, beam)
                fn = apfn1 = os.path.join(
                    self.basedir, "ap_%sb%s.tfs" % (self.name, beam)
                )
                t = getattr(self, "t" + beam).select(n1, n2)
            else:
                t = getattr(self, "t" + beam).select(n1, n2)
            t.plot("x y")
            ffn = "cross_ir%sb%s.png" % (irn, beam)
            pl.title("ir%sb%s" % (irn, beam))
            pl.savefig(ffn)
        return self
