import re

import matplotlib.pyplot as pl

from numpy import *
import numpy as np

from .pydataobj import dataobj
from . import tfsdata

from .poly_fit import poly_fit, poly_print, poly_val

import scipy.interpolate


class StrTable(dataobj):
    scale = 23348.89927

    @classmethod
    def open(cls, fn):
        obj = cls(tfsdata.open(fn))
        return obj

    def get_vars(self, reg):
        rxp = re.compile(reg)
        return sorted(l for l in list(self.keys()) if rxp.match(l))

    def get_kq(self, n):
        return self.get_vars(r"kq[xt]?l?%da?\." % n)

    def get_phases(self):
        out = self.get_vars(r"mu[xy]ip[1-8]b[12]$")
        out += self.get_vars(r"mu[xy]ip[1-8]b[12]_l")
        out += self.get_vars(r"mu[xy]ip[1-8]b[12]_r")
        return out

    def get_betas(self):
        out = self.get_vars(r"bet[xy]ip[1-8]b[12]")
        return out

    def get_disps(self):
        out = self.get_vars(r"d[xy]ip[1-8]b[12]")
        out += self.get_vars(r"dp[xy]ip[1-8]b[12]")
        return out

    def get_acb(self, n, knob="on_sep"):
        out = []
        for n in self.get_vars("acb.*%d\.[lr][1-8].*%s" % (n, knob)):
            if sum(abs(self[n])) > 0:
                out.append(n)
        return sorted(out)

    def get_arcs(self):
        out = self.get_vars(r"kq[fd].a[1-8][1-8]")
        out += self.get_vars(r"kqt[fd].a[1-8][1-8]b[12]")
        out += self.get_vars(r"ks[fd][12].a[1-8][1-8]b[12]")
        out += self.get_vars(r"q[xy]b[12]")
        out += self.get_vars(r"qp[xy]b[12]")
        out += self.get_vars(r"mu[xy][1-8][1-8]b[12]$")
        return sorted(out)

    def plot_arcs(self, n1=None, n2=None, x=None):
        fig = pl.figure("squeeze", figsize=(16, 10))
        fig.clf()
        pl.subplot(4, 4, 1)
        self.plot_mq("kqt?[fd].a81", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 2)
        self.plot_mq("kqt?[fd].a12", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 3)
        self.plot_mq("kqt?[fd].a45", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 4)
        self.plot_mq("kqt?[fd].a56", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 5)
        self.plot_mq("kqt?[fd].a23", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 6)
        self.plot_mq("kqt?[fd].a34", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 7)
        self.plot_mq("kqt?[fd].a67", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 8)
        self.plot_mq("kqt?[fd].a78", n1=n1, n2=n2, x=x)
        pl.subplot(4, 4, 8 + 1)
        self.plot_ms("ks?[fd]..a81", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 2)
        self.plot_ms("ks?[fd]..a12", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 3)
        self.plot_ms("ks?[fd]..a45", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 4)
        self.plot_ms("ks?[fd]..a56", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 5)
        self.plot_ms("ks?[fd]..a23", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 6)
        self.plot_ms("ks?[fd]..a34", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 7)
        self.plot_ms("ks?[fd]..a67", n1=n1, n2=n2, x=x, brho=None)
        pl.subplot(4, 4, 8 + 8)
        self.plot_ms("ks?[fd]..a78", n1=n1, n2=n2, x=x, brho=None)
        pl.tight_layout()
        return self

    def plot_mq(self, vv, n1=None, n2=None, x=None, brho=True):
        if brho is None:
            scale = 1
        else:
            scale = self.scale
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        for k in self.get_vars(vv):
            pl.plot(xv, self[k][n1:n2] * scale, label=k)
            pl.legend(loc=0, frameon=False)
        if brho is None:
            pl.ylabel("k [m${}^{-3}$]")
        else:
            pl.ylabel("k [T/m]")
        return self

    def plot_ms(self, vv, n1=None, n2=None, x=None, brho=True):
        if brho is None:
            scale = 1
        else:
            scale = self.scale
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        for k in self.get_vars(vv):
            pl.plot(xv, self[k][n1:n2] * scale, label=k)
            pl.legend(loc=0, frameon=False)
        if brho is None:
            pl.ylabel("k [m${}^{-3}$]")
        else:
            pl.ylabel("k [T/m${}^2$]")
        return self

    def plot_acb(self, n, knob, n1, n2, x=None, scale=1, brho=None):
        if brho is None:
            scale *= 1e6
        else:
            scale *= brho
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        ks = self.get_acb(n, knob)
        for k in ks:
            pl.plot(xv, self[k][n1:n2] * scale, label=k)
        pl.legend(loc=0, frameon=False)
        if x is not None:
            pl.xlabel(x)
        if brho is None:
            pl.ylabel(r"angle [$\mu$rad]")
        else:
            pl.ylabel("k0l [Tm]")

    def get_triplet(self, trim=True):
        tp = self.get_vars("kqx[123]?[ab]?\.[rl][2815]")
        if trim:
            tp += self.get_vars("ktqx[123]\.[rl][2815]")
        return tp

    def plot_triplet(self, n1, n2, x=None):
        scale = self.scale
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        ks = self.get_triplet()
        for k in ks:
            pl.plot(xv, abs(self[k][n1:n2] * scale), label=k)
        pl.legend(loc=0, frameon=False)
        if x is not None:
            pl.xlabel(x)
        pl.ylabel("k [T/m]")
        a, b = pl.xticks()
        pl.xticks(a[::2])

    def plot_2in1(self, kq, n1, n2, x=None, sign=False, ylab="k [T/m]"):
        scale = self.scale
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        for k in self.get_kq(kq):
            kv = self[k][n1:n2] * scale
            if sign or kv[0] > 0:
                pl.plot(xv, kv, label=k)
            else:
                pl.plot(xv, -kv, label="-" + k)
        pl.legend(loc=0, frameon=False)
        if x is not None:
            pl.xlabel(x)
        pl.ylabel(ylab)
        a, b = pl.xticks()
        pl.xticks(a[::2])

    def plot_ipbeta(self, n1=None, n2=None, x=None):
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        for k in self.get_vars("bet"):
            kv = self[k][n1:n2]
            pl.plot(xv, kv, label=k)
        pl.legend(loc=0, frameon=False)
        if x is not None:
            pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")

    def plot_phase(self, n1=None, n2=None, x=None):
        if x is None:
            xv = arange(len(self[list(self.keys())[0]][n1:n2]))
        else:
            xv = self[x][n1:n2]
        colors = "bgrc"
        for k in self.get_phases():
            kv = self[k][n1:n2]
            pl.plot(xv, kv, color=colors[0], label=k)
            colors = colors[1:] + colors[0]
        pl.legend(loc=0, frameon=False)
        if x is not None:
            pl.xlabel(x)
        pl.ylabel(r"mu [2$\pi$]")
        a, b = pl.xticks()
        pl.xticks(a[::2])

    def plot_squeeze(self, n1=0, n2=None, x=None, plot_phase=True):
        fig = pl.figure("squeeze", figsize=(16, 10))
        fig.canvas.mpl_connect("button_release_event", self.button_press)
        pl.clf()
        if len(self.get_vars("kqx")) > 0:
            pl.subplot(3, 4, 1)
            self.plot_triplet(n1, n2, x=x)
        for n in range(4, 11):
            if len(self.get_kq(n)) > 0:
                pl.subplot(3, 4, n - 2)
                self.plot_2in1(n, n1, n2, x=x, sign=False)
        for n in range(11, 14):
            if len(self.get_kq(n)) > 0:
                pl.subplot(3, 4, n - 2)
                self.plot_2in1(n, n1, n2, x=x, sign=True)
        if plot_phase:
            pl.subplot(3, 4, 12)
            # self.plot_ipbeta(n1,n2,x=x)
            self.plot_phase(n1, n2, x=x)
        # pl.tight_layout()
        self.xvar = x
        return self

    def plot_ir6(self, figname=None, x="scyir5"):
        xv = self[x]
        if figname is None:
            fig = pl.figure(self.filename, figsize=(16, 10))
        else:
            fig = pl.figure(figname, figsize=(16, 10))
        pl.clf()
        pl.subplot(3, 4, 1)
        w = 360
        pl.plot(xv, w * self.refdmuxkickb1_tcdqa, label="TCDQ.A B1")
        pl.plot(xv, w * self.refdmuxkickb1_tcdqb, label="TCDQ.B B1")
        pl.plot(xv, w * self.refdmuxkickb1_tcdqc, label="TCDQ.C B1")
        pl.plot(xv, w * self.refdmuxkickb2_tcdqa, label="TCDQ.A B2")
        pl.plot(xv, w * self.refdmuxkickb2_tcdqb, label="TCDQ.B B2")
        pl.plot(xv, w * self.refdmuxkickb2_tcdqc, label="TCDQ.C B2")
        pl.xlabel(x)
        pl.ylabel(r"$\Delta \mu_x$ MKD-TCDQ [degree]")
        pl.axhline(94, color="k")
        pl.axhline(86, color="k")
        pl.legend()
        pl.subplot(3, 4, 2)
        pl.plot(xv, self.refbxdumpb1, label=r"$\beta_x$ TDE B1")
        pl.plot(xv, self.refbxdumpb2, label=r"$\beta_x$ TDE B2")
        pl.axhline(4000, color="k", label="lim1")
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 3)
        pl.plot(xv, self.refbydumpb1, label=r"$\beta_y$ TDE B1")
        pl.plot(xv, self.refbydumpb2, label=r"$\beta_y$ TDE B2")
        pl.axhline(3200, color="k", label="lim1")
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 4)
        pl.plot(xv, self.refbdumpb1, label=r"$\sqrt{\beta_x\beta_y}$ TDE B1")
        pl.plot(xv, self.refbdumpb2, label=r"$\sqrt{\beta_x\beta_y}$ TDE B2")
        pl.axhline(4500, color="k", label="lim1")
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 5)
        pl.plot(xv, self.refbetxtcdqb1, label=r"$\beta_x$ TCDQ.A B1")
        pl.plot(xv, self.refbetxtcdqb2, label=r"$\beta_x$ TCDQ.A B2")
        # pl.axhline(500,color='k',label='lim')
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 6)
        pl.plot(xv, self.refbetytcdqb1, label=r"$\beta_y$ TCDQ.A B1")
        pl.plot(xv, self.refbetytcdqb2, label=r"$\beta_y$ TCDQ.A B2")
        pl.axhline(145, color="k", label="lim")
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 7)
        pl.plot(xv, self.refbetxtcdsb1, label=r"$\beta_x$ TCDSA.4 B1")
        pl.plot(xv, self.refbetxtcdsb2, label=r"$\beta_x$ TCDSA.4 B2")
        pl.axhline(175, color="k", label="lim at inj.")
        # pl.fill_between([x[0],x[-1]],[200,200],180,color='k',label='lim',alpha=.3)
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 8)
        pl.plot(xv, self.refbetytcdsb1, label=r"$\beta_y$ TCDSA.4 B1")
        pl.plot(xv, self.refbetytcdsb2, label=r"$\beta_y$ TCDSA.4 B2")
        pl.axhline(200, color="k", label="lim")
        pl.xlabel(x)
        pl.ylabel(r"$\beta$ [m]")
        pl.legend()
        pl.subplot(3, 4, 9)
        gaptcdqb1 = 10.1 * sqrt(self.refbetxtcdqb1 * 2.5e-6 / 7000 * 0.938)
        gaptcdqb2 = 10.1 * sqrt(self.refbetxtcdqb2 * 2.5e-6 / 7000 * 0.938)
        pl.plot(xv, 1e3 * gaptcdqb1, label=r"gap TCDQ.4 B1 $10.1\sigma$")
        pl.plot(xv, 1e3 * gaptcdqb2, label=r"gap TCDQ.4 B2 $10.1\sigma$")
        pl.axhline(3, color="k", label="limit 3 mm")
        pl.axhline(4, color="k", label="margin 1 mm")
        pl.xlabel(x)
        pl.ylabel(r"gap [mm] at 10.1$\sigma$")
        pl.legend()
        pl.subplot(3, 4, 10)
        pl.plot(xv, self.refdxtcdqb1, label=r"$D_x$ TCDQ.4 B1")
        pl.plot(xv, self.refdxtcdqb2, label=r"$D_x$ TCDQ.4 B2")
        # pl.plot(xv,self.refdxq4r6b1, label=r'$D_x$ Q4 B1')
        # pl.plot(xv,self.refdxq4r6b2, label=r'$D_x$ Q4 B2')
        pl.xlabel(x)
        pl.ylabel(r"$D_x$ [m]")
        pl.legend()
        pl.subplot(3, 4, 11)
        if not hasattr(self, "ddtct1b1o"):
            self.ddtct1b1o = self.ddtct1b1
            self.ddtct1b2o = self.ddtct1b2
            self.ddtct5b1o = self.ddtct5b1
            self.ddtct5b2o = self.ddtct5b2
        pl.plot(xv, self.ddtct1b1o, color="C1", label=r"MKD-TCT1 B1")
        pl.plot(xv, self.ddtct1b2o, color="C2", label=r"MKD-TCT1 B2")
        pl.plot(xv, self.ddtct5b1o, color="C3", label=r"MKD-TCT5 B1")
        pl.plot(xv, self.ddtct5b2o, color="C4", label=r"MKD-TCT5 B2")
        pl.plot(xv, self.ddtct1b1a, color="C1")
        pl.plot(xv, self.ddtct1b2a, color="C2")
        pl.plot(xv, self.ddtct5b1a, color="C3")
        pl.plot(xv, self.ddtct5b2a, color="C4")
        pl.fill_between(xv, self.ddtct1b1o, self.ddtct1b1a, color="C1", alpha=0.3)
        pl.fill_between(xv, self.ddtct1b2o, self.ddtct1b2a, color="C2", alpha=0.3)
        pl.fill_between(xv, self.ddtct5b1o, self.ddtct5b1a, color="C3", alpha=0.3)
        pl.fill_between(xv, self.ddtct5b2o, self.ddtct5b2a, color="C4", alpha=0.3)
        pl.xlabel(x)
        pl.ylabel(r"$\Delta \mu_x$ [degree]")
        pl.legend()
        pl.subplot(3, 4, 12)
        self["kq5.l6b1"] *= 1.01
        self["kq5.l6b2"] *= 1.01
        self["kq5.r6b1"] *= 1.01
        self["kq5.r6b2"] *= 1.01
        self.plot_2in1(5, 0, None, x="scxir5")
        self["kq5.l6b1"] /= 1.01
        self["kq5.l6b2"] /= 1.01
        self["kq5.r6b1"] /= 1.01
        self["kq5.r6b2"] /= 1.01
        pl.axhline(173, color="k", label="3900 A")
        pl.axhline(175, color="k", linestyle="--", label="3950 A")
        pl.xlabel(x)
        pl.ylabel(r"gradient at 7 TeV [T/m]")
        pl.legend()
        pl.tight_layout()
        return self

    def plot_betsqueeze(
        self, n1=0, n2=None, figname=None, plot_phase=True, semilog=True
    ):
        x = self.get_vars("betxip")[0]
        xv = self[x]
        if figname is None:
            fig = pl.figure(x, figsize=(16, 10))
        else:
            fig = pl.figure(figname, figsize=(16, 10))
        fig.canvas.mpl_connect("button_release_event", self.button_press)
        pl.clf()
        pl.subplot(3, 4, 1)
        if len(self.get_vars("kqx")) > 0:
            self.plot_triplet(n1, n2, x=x)
            if semilog:
                pl.xlim(xv.min(), xv.max())
                pl.semilogx()
        for n in range(4, 11):
            pl.subplot(3, 4, n - 2)
            self.plot_2in1(n, n1, n2, x=x, sign=False)
            if semilog:
                pl.xlim(xv.min(), xv.max())
                pl.semilogx()
        for n in range(11, 14):
            pl.subplot(3, 4, n - 2)
            self.plot_2in1(n, n1, n2, x=x, sign=True)
            if semilog:
                pl.xlim(xv.min(), xv.max())
                pl.semilogx()
        if plot_phase:
            pl.subplot(3, 4, 12)
            self.plot_phase(n1, n2, x=x)
            if semilog:
                pl.xlim(xv.min(), xv.max())
                pl.semilogx()
        pl.tight_layout()
        self.xvar = x
        return self

    def plot_knobs(self, n1=0, n2=None, figname=None, scales=[1, 1]):
        x = self.get_vars("betxip")[0]
        if figname is None:
            fig = pl.figure("knobs", figsize=(16, 10))
        else:
            fig = pl.figure(figname, figsize=(16, 10))
        fig.canvas.mpl_connect("button_release_event", self.button_press)
        pl.clf()
        for ii, (knob, scale) in enumerate(zip(["on_x", "on_sep"], scales)):
            pl.subplot(2, 3, 1 + ii * 3)
            pl.title("%s=%g" % (knob, scale))
            self.plot_acb(1, knob, n1, n2, x=x, scale=scale)
            self.plot_acb(2, knob, n1, n2, x=x, scale=scale)
            self.plot_acb(3, knob, n1, n2, x=x, scale=scale)
            pl.ylim(-100, 100)
            pl.subplot(2, 3, 2 + ii * 3)
            pl.title("%s=%g" % (knob, scale))
            self.plot_acb(4, knob, n1, n2, x=x, scale=scale)
            pl.ylim(-100, 100)
            pl.subplot(2, 3, 3 + ii * 3)
            pl.title("%s=%g" % (knob, scale))
            self.plot_acb(5, knob, n1, n2, x=x, scale=scale)
            if self.get_acb(6):
                self.plot_acb(6, knob, n1, n2, x=x, scale=scale)
            pl.ylim(-100, 100)
        pl.tight_layout()
        self.xvar = x
        return self

    def linint(t, name, val1, val2, x, par="ttt"):
        n1 = where(t[x] == val1)[0][-1]
        n2 = where(t[x] == val2)[0][-1]
        tmp = "%s:=(%g)*(1-%s)+(%g)*(%s);"
        print(tmp % (name, t[name][n1], par, t[name][n2], par))

    def button_press(self, event):
        self.event = event

    def poly_fit_val(
        self, var, order, n1=None, n2=None, param=None, slope0=[], curvep=[]
    ):
        if param is None:
            param = [k for k in list(self.keys()) if k.startswith("betx")][0]
        scale = self.scale
        x = self[param][n1:n2]
        y = self[var][n1:n2]
        x0 = [x[0], x[-1]]
        y0 = [y[0], y[-1]]
        for nn in curvep:
            x0.append(x[nn])
            y0.append(y[nn])
        xp0 = []
        yp0 = []
        for idx in slope0:
            xp0.append(x[idx])
            yp0.append(0)
        pol = poly_fit(order, x, y, x0, y0, xp0, yp0)
        xv = self[param]
        yv = poly_val(pol, xv)
        return yv

    def poly_fit(self, var, order, n1=None, n2=None, param=None, slope0=[], curvep=[]):
        print(var)
        if param is None:
            param = [k for k in list(self.keys()) if k.startswith("betx")][0]
        scale = self.scale
        x = self[param][n1:n2]
        y = self[var][n1:n2]
        x0 = [x[0], x[-1]]
        y0 = [y[0], y[-1]]
        for nn in curvep:
            x0.append(x[nn])
            y0.append(y[nn])
        xp0 = []
        yp0 = []
        for idx in slope0:
            xp0.append(x[idx])
            yp0.append(0)
        print(x, y, x0, y0, xp0, yp0)
        pol = poly_fit(order, x, y, x0, y0, xp0, yp0)
        out = "%s:=%s;" % (var, poly_print(pol, x=param, power="^"))
        if re.match("kt?qx[0-9]?[ab]?\.", var):
            n = 3
        else:
            n = re.match("[kqtxl]+([0-9]+)\.", var)
            if n is None:
                n = 14
                scale = 1
            else:
                n = int(n.groups()[0])
        pl.subplot(3, 4, n - 2)
        xv = self[param]
        yv = poly_val(pol, xv)
        # self[var+"_fit"]=yv
        print("!", " ".join(["%2d" % i for i in sign(diff(yv))]))
        if n < 11:
            yv = abs(yv)
        pl.plot(xv, yv * scale)
        return out

    def poly_fit_all(self, order, param, fn, n1=None, n2=None, slope0=[], curvep=[]):
        out = []
        lst = []
        lst.extend(self.get_triplet())
        lst.extend([kq for n in range(4, 14) for kq in self.get_kq(n)])
        lst.extend(self.get_phases())
        lst.extend(self.get_betas())
        lst.extend(self.get_disps())
        lst.extend(self.get_arcs())
        for vv in lst:
            if param != vv:
                out.append(
                    self.poly_fit(
                        vv,
                        order,
                        param=param,
                        n1=n1,
                        n2=n2,
                        slope0=slope0,
                        curvep=curvep,
                    )
                )
        open(fn, "w").write("\n".join(out))

    def check_slopes(self):
        for kq in range(4, 14):
            for name in self.get_kq(kq):
                v = self[name]
                slopes = sign(diff(v))
                if sum(abs(diff(slopes))) > 0:
                    print(name, " ".join(["%2d" % i for i in slopes]))

    def set_log(self):
        fig = pl.gcf()
        for ax in fig.axes:
            ax.set_xscale("log")
        pl.draw()
        return self

    def set_xlim(self, a, b):
        fig = pl.gcf()
        for ax in fig.axes:
            ax.set_xlim(a, b)
        pl.draw()
        return self

    def savefig(self):
        self.plot_betsqueeze()
        pl.savefig(self.filename.replace(".tfs", ".png"))
        if len(self.get_acb(1)) > 0:
            self.plot_knobs()
            pl.savefig(self.filename.replace(".tfs", "_knobs.png"))
        return self

    def select(self, n1, n2):
        cls = self.__class__
        data = {}
        for key, val in list(self._data.items()):
            if hasattr(val, "dtype"):
                if n2 > n1:
                    data[key] = val[n1:n2].copy()
                else:
                    data[key] = val[n1:n2:-1].copy()
            elif hasattr(val, "copy"):
                data[key] = val.copy()
            else:
                data[key] = val
        return cls(data)

    def interpolate(self, varname, varvalues):
        cls = self.__class__
        data = {}
        xx = self[varname]
        for key, val in list(self._data.items()):
            if hasattr(val, "dtype"):
                newval = scipy.interpolate.interp1d(
                    xx, val, "cubic", assume_sorted=False
                )
                data[key] = newval(varvalues)
            elif hasattr(val, "copy"):
                data[key] = val.copy()
            else:
                data[key] = key
        return cls(data)

    def mk_function(self, variable):
        vv = self[variable]
        xx = np.linspace(0, 1, len(vv), dtype=float)
        return scipy.interpolate.interp1d(xx, vv, "cubic", assume_sorted=False)

    def merge(self, other):
        if self.length() != other.length():
            aa = self.length()
            bb = other.length()
            raise ValueError(f"Tables do not have same length {aa}!={bb}")
        cls = self.__class__
        data = {}
        for key, val in list(self._data.items()):
            if hasattr(val, "copy"):
                data[key] = val.copy()
            else:
                data[key] = val
        for key, val in list(other._data.items()):
            if hasattr(val, "copy"):
                data[key] = val.copy()
            else:
                data[key] = val
        data["col_names"] = self._col_names + other._col_names
        return cls(data)

    def append(self, other):
        cls = self.__class__
        data = {}
        for key, val in list(self._data.items()):
            if hasattr(val, "shape"):
                data[key] = np.r_[val, other._data[key]]
            else:
                data[key] = val
        return cls(data)

    def select(self, keys):
        cls = self.__class__
        data = {}
        for key, val in list(self._data.items()):
            if key.upper() in self._col_names:
                if key in keys:
                    data[key] = val.copy()
            else:
                if hasattr(val, "copy"):
                    data[key] = val.copy()
                else:
                    data[key] = val
        data["col_names"] = [key.upper() for key in keys]
        return cls(data)

    def trim(self, n1, n2):
        cls = self.__class__
        data = {}
        for key, val in list(self._data.items()):
            if key.upper() in self._col_names:
                data[key] = val[n1:n2]
            else:
                if hasattr(val, "copy"):
                    data[key] = val.copy()
                else:
                    data[key] = val
        return cls(data)

    def resample(self, newvalues, varname=None):
        cls = self.__class__
        data = {}
        if varname is None:
            xvalues = np.linspace(0, 1, self.length())
            reverse = False
        else:
            xvalues = self[varname]
            dxmax = max(np.diff(xvalues))
            dxmin = min(np.diff(xvalues))
            if dxmax > 0 and dxmin > 0:
                reverse = False
            elif dxmax < 0 and dxmin < 0:
                reverse = True
            else:
                raise ValueError(f"variable {varname} not monotonus")
        for key, val in list(self._data.items()):
            if key.upper() in self._col_names:
                if reverse:
                    data[key] = np.interp(newvalues, xvalues[::-1], val[::-1])
                else:
                    data[key] = np.interp(newvalues, xvalues, val)
            else:
                if hasattr(val, "copy"):
                    data[key] = val.copy()
                else:
                    data[key] = val
        return cls(data)

    def length(self):
        return len(self[self._col_names[0].lower()])

    def dump_madx(self, val, param, fname):
        out = []
        xx = self[param]
        for k in self._col_names:
            vv = self[k.lower()]
            out.append(f"{k.lower()}={np.interp(val,xx,vv)};")
        out.append(f"{param}={val};")
        open(fname, "w").write("\n".join(out))

    def add_col(self, name, val):
        if len(val) != self.length():
            aa = self.length()
            bb = len(val)
            raise ValueError(f"Column does not have same table length {aa}!={bb}")
        if not name.upper() in self._col_names:
            self._col_names.append(name.upper())
        self[name] = val
