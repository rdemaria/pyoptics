from __future__ import print_function

import re
import os
import gzip
import time


import matplotlib.pyplot as pl
import matplotlib
from matplotlib.animation import FuncAnimation
import numpy as np
import scipy


from .utils import mystr as _mystr
from .utils import pyname
from collections import namedtuple
from .pydataobj import dataobj
from . import tfsdata
from .survey import rot_mad, get_w_from_angles
from .tablemixin import TableMixIn


def rng(x, a, b):
    "return (x<b) & (x>a)"
    return (x < b) & (x > a)


def errors_getmn(self, errorname="b6"):
    """returns list of (m,n) for multiople error errorname"""
    kind = errorname[0]
    order = int(errorname[1:])
    return getmn(order=order, kind=kind)


def getmn(order, kind="b"):
    """returns list of (m,n) of * resonances of order o
    with * = 't': all resonances
             'a': skew multipoles n=odd
             'b': normal multipoles n=even
             's': sum resonances (m>0,n>0), loss of beam
             'd': difference resonances (m<0,n>0) or (m>0,n<0), exchange between planes
    """

    """Return resonances given the order and type"""
    out = []
    if "t" in kind:
        kind = "ab"
    for m in range(0, order + 1):
        n = order - m
        if "b" in kind and n % 2 == 0 or m == 0:
            out.append((m, n))
            if n > 0:
                out.append((m, -n))
        if "a" in kind and n % 2 == 1 and m > 0:
            out.append((m, n))
            if n > 0:
                out.append((m, -n))
        if "s" in kind and (n > 0 and m > 0):
            out.append((m, n))
        if "d" in kind and (n > 0 and m > 0):
            out.append((m, -n))
    return list(set(out))


class optics(dataobj, TableMixIn):
    _is_s_begin = False
    _name_char = 16
    _entry_char = 12
    _entry_prec = 3

    @classmethod
    def open(cls, fn):
        return cls(tfsdata.open(fn))

    def save(self, fn, floatfmt="%20.9f"):
        tfsdata.save(self._data, fn, floatfmt)

    def __init__(self, data={}, idx=False):
        if hasattr(data, "col_names"):
            self._data["_col_names"] = data.col_names()
        self.update(data)
        if hasattr(data, "summary"):
            self.header = data.summary
        if hasattr(data, "name"):
            self.name = np.array([a.split(":")[0] for a in self.name])
        self._fdate = 0

    def col_names(self):
        if "_col_names" in self._data:
            return self._data["_col_names"]
        else:
            return self._data.col_names()

    def copy(self):
        data = {}
        for k, v in list(self._data.items()):
            if hasattr(v, "copy"):
                vv = v.copy()
            #      elif  hasattr(v,'__getitem__'):
            #        vv=v[:]
            else:
                vv = v
            data[k] = vv
        return optics(data)

    def reload(self):
        if "filename" in self._data:
            fdate = os.stat(self.filename).st_ctime
            if fdate > self._fdate:
                self._data = tfsdata.open(self.filename)
                self._fdate = fdate
                return True
        return False

    def get_idx(self, name=None, count=0):
        if type(name) is str:
            return np.where(self.name == name)[0][count]
        else:
            return count

    def row(self, index):
        return {cc: self[cc][index] for cc in self.col_names()()}

    def twissdata(self, location, data):
        idx = np.where(self.pattern(location))[0][-1]
        out = dict(location=location)
        for name in data.split():
            vec = self.__dict__.get(name)
            if vec is None:
                out[name] = 0
            else:
                out[name] = vec[idx]
        out["sequence"] = self.header.get("sequence")
        return out

    def range(self, pat1, pat2):
        """return a mask relative to range"""
        try:
            id1 = np.where(self.pattern(pat1))[0][-1]
        except IndexError:
            raise ValueError("%s pattern not found in table" % pat1)
        try:
            id2 = np.where(self.pattern(pat2))[0][-1]
        except IndexError:
            raise ValueError("%s pattern not found in table" % pat2)
        out = np.zeros(len(self.name), dtype=bool)
        if id2 > id1:
            out[id1 : id2 + 1] = True
        else:
            out[id1:] = True
            out[: id2 + 1] = True
        return out

    def plot(
        self,
        yl="",
        yr="",
        x="s",
        idx=slice(None),
        clist="k r b g c m",
        lattice=True,
        newfig=True,
        pre=None,
        autoupdate=False,
        ip_label=False,
    ):
        out = qdplot(
            self,
            x=x,
            yl=yl,
            yr=yr,
            idx=idx,
            lattice=lattice,
            newfig=newfig,
            clist=clist,
            pre=pre,
        )
        self._plot = out
        if ip_label:
            self.set_xaxis_ip()
        if autoupdate:
            if "wx" in matplotlib.get_backend().lower():
                self._plot.wx_autoupdate()
            else:
                self._plot.ani_autoupdate()
        #    self._target.append(out)
        return self

    def plotbeta(self, **nargs):
        return self.plot("betx bety", "dx dy", **nargs)

    def plotsigma(self, emit=2.5e-6 / 7000 * 0.938, deltap=1.1e-4, **nargs):
        self.sigx = np.sqrt(self.betx * emit) * 1000
        self.sigy = np.sqrt(self.bety * emit) * 1000
        self.sigdx = self.dx * deltap * 1000
        self.plot("sigx sigy sigdx", **nargs)
        ya, yb = pl.ylim()
        pl.twinx()
        bmax = max(self.betx.max(), self.bety.max())
        rng = list(range(0, int(np.ceil(np.log10(bmax))) + 1))
        bval = np.array([n * 10**dd for dd in rng for n in [1, 2, 5]])
        bval = bval[bval < bmax]
        pl.ylim(ya, yb)
        return self

    def plotcross(self, **nargs):
        return self.plot("x y", "dx dy", **nargs)

    def plottune(self, newfig=True, **nargs):
        q4x, q3x, q2x, q1x, q0x = scipy.polyfit(self.deltap, self.q1, 4)
        q4y, q3y, q2y, q1y, q0y = scipy.polyfit(self.deltap, self.q2, 4)
        qx = (self.q1 - self.q1.round())[abs(self.deltap) < 1e-15][0]
        qy = (self.q2 - self.q2.round())[abs(self.deltap) < 1e-15][0]
        fmt = r"$%s=%4.2f %+4.2f \delta"
        fmt += r"%+4.2f \frac{\delta^2}{2\cdot10^{-3}}"
        fmt += r"%+4.2f \frac{\delta^3}{6\cdot10^{-6}}$"
        fmtx = fmt % ("Q_x", q0x, q1x, q2x * 2e-3, q3x * 6e-6)
        fmty = fmt % ("Q_y", q0y, q1y, q2y * 2e-3, q3y * 6e-6)
        if newfig:
            pl.figure()
        pl.title(r"Tune vs $\delta$")
        pl.xlabel("$\delta$")
        pl.ylabel("Fractional tune")
        pl.plot(self.deltap, self.q1 - self.q1.round(), label=fmtx, **nargs)
        pl.plot(self.deltap, self.q2 - self.q2.round(), label=fmty, **nargs)
        # pl.text(0.0,qx,r"$Q_x$")
        # pl.text(0.0,qy,r"$Q_y$")
        pl.grid(True)
        pl.legend(loc=0)
        return self

    def plotbetabeat(self, t1, dp="0.0003", **nargs):
        pl.title(r"$\rm{Beta beat: 1 - \beta(\delta=%s)/\beta(\delta=0)}$" % dp)
        pl.ylabel(r"$\Delta\beta/\beta$")
        pl.xlabel(r"$s$")
        pl.plot(
            self.s, 1 - t1.betx / self.betx, label=r"$\Delta\beta_x/\beta_x$", **nargs
        )
        pl.plot(
            self.s, 1 - t1.bety / self.bety, label=r"$\Delta\beta_y/\beta_y$", **nargs
        )
        pl.grid(True)
        pl.legend()
        return self

    def plotw(self, lbl="", **nargs):
        pl.title(r"Chromatic function: %s" % lbl)
        # pl.ylabel(r"$w=(\Delta\beta/\beta)/\delta$")
        pl.ylabel(r"$w$")
        pl.xlabel(r"$s$")
        pl.plot(self.s, self.wx, label=r"$w_x$", **nargs)
        pl.plot(self.s, self.wy, label=r"$w_y$", **nargs)
        pl.grid(True)
        pl.legend()
        return self

    def plotap(self, ap=None, nlim=30, ref=12.6, newfig=True, eref=None, **nargs):
        if ap is None:
            apfn = self.filename.replace("twiss", "ap")
            ap = optics.open(apfn)
        if eref is not None:
            ap.s -= ap.s[ap // eref]
            self.s -= self.s[self // eref]
        self.ss = ap.s
        self.n1 = ap.n1
        self = self.plot(x="ss", yl="n1", newfig=newfig, **nargs)
        p = self._plot
        p.figure.gca().plot(self.ss, self.ss * 0 + ref)
        p.figure.gca().set_ylim(0, nlim)
        p.figure.canvas.draw()
        self._plot = p
        self.ap = ap
        return self

    def mk_betamax(self):
        self.betxmax = np.zeros_like(self.betx)
        self.betymax = np.zeros_like(self.bety)
        bufname = ""
        for i in range(len(self.k1l)):
            k1l = self.k1l[i]
            l = self.l[i]
            betx = self.betx[i]
            alfx = self.alfx[i]
            bety = self.bety[i]
            alfy = self.alfy[i]
            if l > 0:  # thick
                alfxm = self.alfx[i - 1]
                alfym = self.alfy[i - 1]
                betxm = self.betx[i - 1]
                betym = self.bety[i - 1]
                if k1l >= 0:  # focus
                    if alfxm < 0 and alfx > 0:
                        self.betxmax[i] = betx + alfx**2 / betx / k1l * l
                    else:
                        self.betxmax[i] = max(betxm, betx)
                    self.betymax[i] = max(bety, betym)
                else:  # defocu
                    if alfym < 0 and alfy > 0:
                        self.betymax[i] = bety + alfy**2 / bety / (-k1l) * l
                    else:
                        self.betymax[i] = max(bety, betym)
                    self.betxmax[i] = max(betxm, betx)
            elif abs(k1l) > 0:  # thin
                name = self.name[i].split("..")[0]
                if bufname != name:
                    if bufname != "":
                        idx = np.where(self.name == bufname)[0][0]
                        self.betxmax[idx] = max(bufx)
                        self.betymax[idx] = max(bufy)
                    bufx = [self.betx[i]]
                    bufy = [self.bety[i]]
                    bufname = name
                else:
                    bufx.append(self.betx[i])
                    bufy.append(self.bety[i])
        return self

    def maxbety(f):
        return f.bety + f.alfy**2 / f.bety / abs(f.k1l / f.l)

    def maxbety(f):
        return f.bety + f.alfy**2 / f.bety / abs(f.k1l / f.l)

    # def chromx(f):
    #  if not hasattr(f,'k1l'):
    #    f.k1l=f.k1l
    #  return -sum(f.k1l*f.betx)/4/pi
    # def chromy(f):
    #  if not hasattr(f,'k1l'):
    #    f.k1l=f.k1l
    #  return sum(f.k1l*f.bety)/4/pi
    def ndx(self):
        return self.dx / np.sqrt(self.betx)

    def ndpx(self):
        return self.dpx * np.sqrt(self.betx) + self.dx / np.sqrt(self.betx) * self.alfx

    def alphac(self):
        return sum(self("dx*kn0l")) / sum(self.l)

    def gammatr(self):
        af = self._alphac()
        if af > 0:
            return np.sqrt(1 / af)
        else:
            return -np.sqrt(-1 / af)

    def transferMatrix(self, i1=0, i2=-1, plane="x"):
        """Return the transfer matrix from position i1 to position i2
        see Y.Lee 2.68 pag 53 for definition
        """
        B2 = self.normMat(i2, plane=plane)
        B1 = self.normMat(i1, plane=plane)
        psi = 2 * np.pi * (self["mu" + plane][i2] - self["mu" + plane][i1])
        R = np.array([[np.cos(psi), np.sin(psi)], [-np.sin(psi), np.cos(psi)]])
        return np.dot(np.dot(B2, R), np.linang.inv(B1))

    def normMat(self, i, plane="x"):
        beta = self["bet" + plane][i]
        alpha = self["alf" + plane][i]
        return np.array(
            [[np.sqrt(beta), 0], [-alpha / np.sqrt(beta), 1 / np.sqrt(beta)]]
        )

    def mk_intkbeta(self, on_sext=True):
        self.intkbetx = self.k1l * 0.0
        self.intkbety = self.k1l * 0.0
        for i in range(len(self.k1l)):
            kl = self.k1l[i]
            k2l = self.k2l[i]
            l = self.l[i]
            betx = self.betx[i]
            bety = self.bety[i]
            intkbetx = intkbety = 0
            if abs(kl) > 0:
                if l == 0:
                    intkbetx = kl * betx
                    intkbety = kl * bety
                else:
                    alfx = self.alfx[i]
                    alfy = self.alfy[i]
                    gamx = (1 + alfx**2) / betx
                    gamy = (1 + alfy**2) / bety
                    k = kl / l
                    ak = np.abs(k)
                    rk = np.sqrt(ak)
                    rkl = rk * l
                    crkl = np.cos(rkl)
                    srkl = -np.sin(rkl)
                    chrkl = np.cosh(rkl)
                    shrkl = -np.sinh(rkl)
                    if k > 0:  # backtrack
                        r11 = crkl
                        r12 = srkl / rk
                        r21 = -srkl * rk
                        r33 = chrkl
                        r34 = shrkl / rk
                        r43 = shrkl * rk
                    else:
                        r33 = crkl
                        r34 = srkl / rk
                        r43 = -srkl * rk
                        r11 = chrkl
                        r12 = shrkl / rk
                        r21 = shrkl * rk
                    r22 = r11
                    r44 = r33
                    betx0 = r11**2 * betx - 2.0 * r11 * r12 * alfx + r12**2 * gamx
                    bety0 = r33**2 * bety - 2.0 * r33 * r34 * alfy + r34**2 * gamy
                    alfx0 = (
                        -r11 * r21 * betx
                        + (1.0 + 2.0 * r12 * r21) * alfx
                        - r12 * r22 * gamx
                    )
                    alfy0 = (
                        -r33 * r43 * bety
                        + (1.0 + 2.0 * r34 * r43) * alfy
                        - r34 * r44 * gamy
                    )
                    gamx0 = (1.0 + alfx0**2) / betx0
                    gamy0 = (1.0 + alfy0**2) / bety0
                    if k > 0:
                        ax = 0.5 * (l + 0.5 / rk * np.sin(2.0 * rkl))
                        bx = 0.25 / ak * (1.0 - np.cos(2.0 * rkl))
                        cx = 0.5 / ak * (l - 0.5 / rk * np.sin(2.0 * rkl))
                        ay = 0.5 * (l + 0.5 / rk * np.sinh(2.0 * rkl))
                        by = -0.25 / ak * (1.0 - np.cosh(2.0 * rkl))
                        cy = -0.5 / ak * (l - 0.5 / rk * np.sinh(2.0 * rkl))
                    else:
                        ay = 0.5 * (l + 0.5 / rk * np.sin(2.0 * rkl))
                        by = 0.25 / ak * (1.0 - np.cos(2.0 * rkl))
                        cy = 0.5 / ak * (l - 0.5 / rk * np.sin(2.0 * rkl))
                        ax = 0.5 * (l + 0.5 / rk * np.sinh(2.0 * rkl))
                        bx = -0.25 / ak * (1.0 - np.cosh(2.0 * rkl))
                        cx = -0.5 / ak * (l - 0.5 / rk * np.sinh(2.0 * rkl))
                    # if (self//'MQX.*L5')[i]:
                    #  print kl,l,k,ak,rk,rkl
                    #  print betx ,bety ,alfx ,alfy ,gamx ,gamy
                    #  print r11,r12,r21,r22,r33,r34,r43,r44
                    #  print betx0,bety0,alfx0,alfy0,gamx0,gamy0
                    #  print ax,bx,cx,ay,by,cy
                    intkbetx = k * (ax * betx0 - 2.0 * bx * alfx0 + cx * gamx0)
                    intkbety = -k * (ay * bety0 - 2.0 * by * alfy0 + cy * gamy0)
            elif abs(k2l) > 0:
                dx = self.dx[i]
                dy = self.dy[i]
                intkbetx = -k2l * dx * betx
                intkbety = k2l * dx * bety
            self.intkbetx[i] = intkbetx
            self.intkbety[i] = intkbety
        return self

    def mk_AB(self):
        """Values
        qpp1: (k2 D -k1)' beta
        qpp2: k2 D''
        qpp3: k1 beta'
        qpp4: k2 D beta'
        """
        betx, bety, dx = self.betx, self.bety, self.dx
        self.Bx = betx * self.wx * np.cos(2 * np.pi * self.phix)
        self.Ax = betx * self.wx * np.sin(2 * np.pi * self.phix)
        self.By = bety * self.wy * np.cos(2 * np.pi * self.phiy)
        self.Ay = bety * self.wy * np.sin(2 * np.pi * self.phiy)
        k1l = self.k1l
        k2l = self.k2l
        self.qp1x = 1.0 / 4 / np.pi * np.cumsum(-betx * k1l)
        self.qp1y = 1.0 / 4 / np.pi * np.cumsum(bety * k1l)
        self.qp2x = 1.0 / 4 / np.pi * np.cumsum(betx * k2l * dx)
        self.qp2y = 1.0 / 4 / np.pi * np.cumsum(-bety * k2l * dx)
        self.qpx = self.qp1x + self.qp2x
        self.qpy = self.qp1y + self.qp2y
        self.qpp1x = -2 * self.qpx
        self.qpp1y = -2 * self.qpy
        self.qpp2x = 1.0 / 2 / np.pi * np.cumsum(k2l * self.ddx * betx)
        self.qpp2y = 1.0 / 2 / np.pi * np.cumsum(-k2l * self.ddx * bety)
        self.qpp3x = 1.0 / 4 / np.pi * np.cumsum(-k1l * self.Bx)
        self.qpp3y = 1.0 / 4 / np.pi * np.cumsum(k1l * self.By)
        self.qpp4x = 1.0 / 4 / np.pi * np.cumsum(k2l * dx * self.Bx)
        self.qpp4y = 1.0 / 4 / np.pi * np.cumsum(-k2l * dx * self.By)
        self.qppx = self.qpp1x + self.qpp2x + self.qpp3x + self.qpp4x
        self.qppy = self.qpp1y + self.qpp2y + self.qpp3y + self.qpp4y
        qp1x = self.qp1x[-1]
        qp2x = self.qp2x[-1]
        qpx = self.qpx[-1]
        qpp1x = self.qpp1x[-1]
        qpp2x = self.qpp2x[-1]
        qpp3x = self.qpp3x[-1]
        qpp4x = self.qpp4x[-1]
        qppx = self.qppx[-1]
        qp1y = self.qp1y[-1]
        qp2y = self.qp2y[-1]
        qpy = self.qpy[-1]
        qpp1y = self.qpp1y[-1]
        qpp2y = self.qpp2y[-1]
        qpp3y = self.qpp3y[-1]
        qpp4y = self.qpp4y[-1]
        qppy = self.qppy[-1]
        print("Qx' = %10g%+10g = %10g" % (qp1x, qp2x, qpx))
        print("Qy' = %10g%+10g = %10g" % (qp1y, qp2y, qpy))
        print("Qx''= %10g%+10g%+10g%+10g = %10g" % (qpp1x, qpp2x, qpp3x, qpp4x, qppx))
        print("Qy''= %10g%+10g%+10g%+10g = %10g" % (qpp1y, qpp2y, qpp3y, qpp4y, qppy))
        self.qp1x = qp1x
        self.qp2x = qp2x
        self.qpx = qpx
        self.qpp1x = qpp1x
        self.qpp2x = qpp2x
        self.qpp3x = qpp3x
        self.qpp4x = qpp4x
        self.qppx = qppx
        self.qp1y = qp1y
        self.qp2y = qp2y
        self.qpy = qpy
        self.qpp1y = qpp1y
        self.qpp2y = qpp2y
        self.qpp3y = qpp3y
        self.qpp4y = qpp4y
        self.qppy = qppy
        return self

    def interp(self, snew, namenew=None, sname="s"):
        "Interpolate with piecewise linear all columns using a new s coordinate"
        for cname in self.col_names():
            if cname != sname and np.isreal(self[cname][0]):
                self[cname] = np.interp(snew, self[sname], self[cname])
        self[sname] = snew
        self.name = namenew

    def _first_idx(self, name):
        if type(name) is str:
            name = np.where(self.name == name)[0][0]
        return name

    def _iter_columns(self):
        ln = len(self.name)
        for k, v in list(self._data.items()):
            if hasattr(v, "__len__") and len(v) == ln:
                yield k, v

    def cycle(self, name, reorder=True):
        idx = self._first_idx(name)
        for vn in ["s", "mux", "muy", "phix", "phiy"]:
            if vn in self:
                v = self[vn]
                vm = v[-1]
                v -= v[idx]
                if reorder:
                    v[:idx] += vm
        if reorder:
            for vn, v in self._iter_columns():
                v = self[vn]
                self[vn] = np.concatenate([v[idx:], v[:idx]])
        if hasattr(self, "ap"):
            self.ap.cycle(idx, reorder=reorder)
        return self

    def center(self, ref):
        idx = np.where(self // ref)[0][0]
        if self.header["type"] in ["TWISS", "APERTURE"]:
            for vn in ["s", "mux", "muy", "phix", "phiy"]:
                if vn in self:
                    v = self[vn]
                    v -= v[idx]
        elif self.header["type"] == "SURVEY":
            theta0 = self.theta[idx]
            c0 = np.cos(theta0)
            s0 = np.sin(theta0)
            x0 = self.x[idx]
            y0 = self.y[idx]
            z0 = self.z[idx]
            xx = self.x - x0
            yy = self.y - y0
            zz = self.z - z0
            xxx = xx * c0 - zz * s0
            zzz = xx * s0 + zz * c0
            self.x = xxx
            self.z = zzz
            self.s -= self.s[idx]
        return self

    def select(self, a, b, shift=True):
        a = self._first_idx(a)
        b = self._first_idx(b)
        data = {}
        ln = len(self.name)
        for k, v in list(self._data.items()):
            if hasattr(v, "__len__") and len(v) == ln:
                vv = v[a : b + 1]
            elif hasattr(v, "copy"):
                vv = v.copy()
            #      elif hasattr(v,'__getitem__'):
            #        vv=v[:]
            else:
                vv = v
            data[k] = vv
        if shift:
            for vn in ["s", "mux", "muy", "phix", "phiy"]:
                if vn in data:
                    data[vn] -= data[vn][0]
        return optics(data)

    def append(self, t):
        data = {}
        for k, v in list(self._data.items()):
            if k in self.col_names():
                data[k] = np.concatenate([v, t[k]])
            else:
                data[k] = v
        return optics(data)

    def resize(self, nn):
        data = {}
        for k, v in list(self._data.items()):
            if k in self.col_names():
                data[k] = np.zeros(nn, dtype=v.dtype)
            else:
                data[k] = v
        return optics(data)

    def errors_add(self, error_table):
        """Add error columns"""
        klist = []
        for k, val in list(error_table.items()):
            if k.startswith("k") and sum(abs(val)) > 0:
                klist.append([k, val])
                self[k] = self.get(k, np.zeros(len(self.name)))
        for idxerror, name in enumerate(error_table["name"]):
            idxself = np.where(self.name == name)[0]
            for k, val in klist:
                self[k][idxself] += val[idxerror]
        return self

    def drvterm(t, m=0, n=0, p=0, q=0):
        dv = t.betx ** (abs(m) / 2.0) * t.bety ** (abs(n) / 2.0)
        dv = dv * np.exp(+2j * np.pi * ((m - 2 * p) * t.mux + (n - 2 * q) * t.muy))
        return dv

    def errors_kvector(self, i, maxorder=10):
        rng = list(range(maxorder))
        kn, ks = [], []
        for n in rng:
            kname = "k%dl" % n
            if kname in self:
                kn.append(self[kname][i])
            else:
                kn.append(0.0)
            kname = "k%dsl" % n
            if kname in self:
                ks.append(self[kname][i])
            else:
                ks.append(0.0)
        return kn, ks

    def errors_ktob(self, maxorder=6):
        nelem = len(self.name)
        rng = list(range(maxorder))
        xx = self.x
        yy = self.y
        for n in rng:
            self["b%d" % (n + 1)] = np.zeros(nelem)
            self["a%d" % (n + 1)] = np.zeros(nelem)
        for i in range(nelem):
            kn, ks = self.errors_kvector(i, maxorder)
            x = xx[i]
            y = yy[i]
            cn = k2b(kn, ks, x, y)
            for ib, b in enumerate(cn):
                self["b%d" % (ib + 1)][i] = b.real
                self["a%d" % (ib + 1)][i] = b.imag
        return self

    def errors_detuning(self, ex, ey, xs, ys, order):
        nelem = len(self.name)
        Jx = ex / 2 * xs**2
        Jy = ey / 2 * ys**2
        Dq = [Dq2, Dq4, Dq6, Dq8, Dq10, Dq12, Dq14, Dq16][order / 2 - 1]
        betx = self.betx
        bety = self.bety
        bnn = self["b%d" % order]
        dqxx = np.zeros(nelem)
        dqyy = np.zeros(nelem)
        for i in range(nelem):
            Bx = betx[i]
            By = bety[i]
            bn = bnn[i]
            dqx, dqy = Dq(bn, Bx, Jx, By, Jy)
            dqxx[i] = dqx
            dqyy[i] = dqy
            # print bn,Bx,Jx,By,Jy,dqx,dqx
        self["DQx%d" % order] = dqxx
        self["DQy%d" % order] = dqyy
        return self

    def errors_footprint(
        self,
        ex=3.75e-6 / 450 * 0.938,
        ey=3.75e-6 / 450 * 0.938,
        nsigma=12,
        nangles=7,
        orders=[4, 6],
        wp=(0.28, 0.31),
        label="footprint",
    ):
        x, y = mk_grid(nsigma, nangles)
        tunx = []
        tuny = []
        for order in orders:
            if order % 2 == 1 or order < 2 or order > 8:
                print("Order supported are 2,4,6,8,10,12,14,16")
        for xs, ys in zip(x, y):
            dqx = wp[0]
            dqy = wp[1]
            for order in orders:
                self.errors_detuning(ex, ey, xs, ys, order)
                dqx += sum(self["DQx%d" % order])
                dqy += sum(self["DQy%d" % order])
            tunx.append(dqx)
            tuny.append(dqy)
        return Footprint(x, y, tunx, tuny, nsigma, nangles)

    def set_xaxis_ip(self):
        idx = self // "IP.$"
        sl = self.s[idx]
        ns = self.name[idx]
        ax = pl.gca()
        ax.set_xticks(sl)
        ax.set_xticklabels(ns)

    def idx_from_namelist(self, namelist):
        iilist = 0
        currname = namelist[iilist]
        out = []
        for ii, name in enumerate(self.name):
            if name == currname:
                iilist += 1
                if iilist < len(namelist):
                    currname = namelist[iilist]
                out.append(ii)
        return out

    def cox(self, elem):
        el = np.where(self // elem)[0][0]
        mu0 = self.mux[el]
        bet0 = self.betx[el]
        pq0 = np.pi * self.header["q1"]
        return (
            0.5
            / np.sin(pq0)
            * np.sqrt(bet0 * self.betx)
            * np.cos(2 * np.pi * abs(mu0 - self.mux) - pq0)
        )

    def coy(self, elem):
        el = np.where(self // elem)[0][0]
        mu0 = self.muy[el]
        bet0 = self.bety[el]
        pq0 = np.pi * self.header["q2"]
        return (
            0.5
            / np.sin(pq0)
            * np.sqrt(bet0 * self.bety)
            * np.cos(2 * np.pi * abs(mu0 - self.muy) - pq0)
        )

    def get_rotmat(self, i):
        return rot_mad(self.theta[i], self.phi[i], self.psi[i])

    def get_pos(self, i):
        return np.array([self.x[i], self.y[i], self.z[i]])


def _mylbl(d, x):
    return d.get(x, r"$%s$" % x)


class qdplot(object):
    lglabel = {
        "betx": r"$\beta_x$",
        "bety": r"$\beta_y$",
        "dx": r"$D_x$",
        "dy": r"$D_y$",
        "mux": r"$\mu_x$",
        "muy": r"$\mu_y$",
        "Ax": "$A_x$",
        "Ay": "$A_y$",
        "Bx": "$B_x$",
        "By": "$B_y$",
        "wx": "$w_x$",
        "wy": "$w_y$",
        "sigx": r"$\sigma_x=\sqrt{\beta_x \epsilon}$",
        "sigy": r"$\sigma_y=\sqrt{\beta_y \epsilon}$",
        "sigdx": r"$\sigma_{D_x}=D_x \delta$",
        "n1": r"Aperture [$\sigma$]",
    }

    axlabel = {
        "s": r"$s [m]$",
        "ss": r"$s [m]$",
        "betx": r"$\beta [m]$",
        "bety": r"$\beta [m]$",
        "mux": r"$\mu/(2 \pi)$",
        "muy": r"$\mu/(2 \pi)$",
        "dx": r"$D [m]$",
        "dy": r"$D [m]$",
        "x": r"$co [m]$",
        "y": r"$co [m]$",
        "sigx": r"$\sigma$ [mm]",
        "sigy": r"$\sigma$ [mm]",
        "sigdx": r"$\sigma$ [mm]",
        "n1": r"Aperture [$\sigma$]",
    }
    autoupdate = []

    def ani_autoupdate(self):
        from matplotlib.animation import FuncAnimation

        self._ani = FuncAnimation(self.figure, self.update, blit=False, interval=1000)

    def ani_stopupdate(self):
        del self._ani

    @classmethod
    def on_updated(cls, fun):
        cls.on_update = fun

    def __init__(
        self,
        t,
        x="",
        yl="",
        yr="",
        idx=slice(None),
        clist="k r b g c m",
        lattice=None,
        newfig=True,
        pre=None,
    ):
        yl, yr, clist = list(map(str.split, (yl, yr, clist)))
        #    timeit('Init',True)
        self.color = {}
        self.left = None
        self.right = None
        self.lattice = None
        self.pre = None
        self.t, self.x, self.yl, self.yr, self.idx, self.clist = (
            t,
            x,
            yl,
            yr,
            idx,
            clist,
        )
        for i in self.yl + self.yr:
            self.color[i] = self.clist.pop(0)
            self.clist.append(self.color[i])
        if newfig is True:
            self.figure = pl.figure()
        elif newfig is False:
            self.figure = pl.gcf()
            self.figure.clf()
        else:
            self.figure = newfig
            self.figure.clf()
        if lattice:
            self.lattice = self._new_axes()
            #      self.lattice.set_autoscale_on(False)
            self.lattice.yaxis.set_visible(False)
        if yl:
            self.left = self._new_axes()
            #      self.left.set_autoscale_on(False)
        if yr:
            self.right = self._new_axes()
            #      self.right.set_autoscale_on(False)
            self.left.yaxis.set_label_position("right")
            self.left.yaxis.set_ticks_position("right")

        #    timeit('Setup')
        self.run()
        if lattice is not None:
            self.lattice.set_autoscale_on(False)
        if yl:
            self.left.set_autoscale_on(False)
            self.left.yaxis.set_label_position("left")
            self.left.yaxis.set_ticks_position("left")
        if yr:
            self.right.set_autoscale_on(False)
            self.right.yaxis.set_label_position("right")
            self.right.yaxis.set_ticks_position("right")

    #    timeit('Update')
    def _new_axes(self):
        if self.figure.axes:
            ax = self.figure.axes[-1]
            out = self.figure.add_axes(ax.get_position(), sharex=ax, frameon=False)
        else:
            # adjust plot dimensions
            out = self.figure.add_axes([0.17, 0.12, 0.6, 0.8])
        return out

    def __repr__(self):
        return object.__repr__(self)

    def _trig(self):
        print("optics trig")
        self.run()

    def update(self, *args):
        if hasattr(self.t, "reload"):
            if self.t.reload():
                self.run()
                return self
        return False

    #  def _wx_callback(self,*args):
    #    self.update()
    #    wx.WakeUpIdle()
    #
    #  def autoupdate(self):
    #    if pl.rcParams['backend']=='WXAgg':
    #      wx.EVT_IDLE.Bind(wx.GetApp(),wx.ID_ANY,wx.ID_ANY,self._wx_callback)
    #    return self
    #
    #  def stop_update(self):
    #    if pl.rcParams['backend']=='WXAgg':
    #      wx.EVT_IDLE.Unbind(wx.GetApp(),wx.ID_ANY,wx.ID_ANY,self._callback)
    #
    #  def __del__(self):
    #    if hasattr(self,'_callback'):
    #      self.stop_update()

    def run(self):
        #    print 'optics run'
        self.ont = self.t
        self.xaxis = getattr(self.ont, self.x)[self.idx]
        is_ion = pl.isinteractive()
        pl.interactive(False)
        self.lines = []
        self.legends = []
        #    self.figure.lines=[]
        #    self.figure.patches=[]
        #    self.figure.texts=[]
        #    self.figure.images = []
        self.figure.legends = []

        if self.lattice:
            self.lattice.clear()
            self._lattice(["k0l", "kn0l", "angle"], "#a0ffa0", "Bend h")
            self._lattice(["ks0l"], "#ffa0a0", "Bend v")
            self._lattice(["kn1l", "k1l"], "#a0a0ff", "Quad")
            self._lattice(["hkick"], "#e0a0e0", "Kick h")
            self._lattice(["vkick"], "#a0e0e0", "Kick v")
            self._lattice(["kn2l", "k2l"], "#e0e0a0", "Sext")
        if self.left:
            self.left.clear()
            for i in self.yl:
                self._column(i, self.left, self.color[i])
        if self.right:
            self.right.clear()
            for i in self.yr:
                self._column(i, self.right, self.color[i])
        ca = self.figure.gca()
        ca.set_xlabel(_mylbl(self.axlabel, self.x))
        ca.set_xlim(min(self.xaxis), max(self.xaxis))
        self.figure.legend(self.lines, self.legends, loc="upper right")
        ca.grid(True)
        #    self.figure.canvas.mpl_connect('button_release_event',self.button_press)
        self.figure.canvas.mpl_connect("pick_event", self.pick)
        pl.interactive(is_ion)
        self.figure.canvas.draw()
        if hasattr(self, "on_run"):
            self.on_run(self)

    def pick(self, event):
        pos = np.array([event.mouseevent.x, event.mouseevent.y])
        name = event.artist.elemname
        prop = event.artist.elemprop
        value = event.artist.elemvalue
        print("\n %s.%s=%s" % (name, prop, value), end=" ")

    #  def button_press(self,mouseevent):
    #    rel=np.array([mouseevent.x,mouseevent.y])
    #    dx,dy=self.pickpos/rel
    #    print 'release'
    #    self.t[self.pickname][self.pickprop]*=dy
    #    self.t.track()
    #    self.update()

    def _lattice(self, names, color, lbl):
        #    timeit('start lattice %s' % names,1)
        vd = 0
        sp = self.lattice
        s = self.ont.s
        l = self.ont.l
        for i in names:
            myvd = self.ont.__dict__.get(i, None)
            if myvd is not None:
                vdname = i
                vd = myvd[self.idx] + vd
        if np.any(vd != 0):
            m = np.abs(vd).max()
            if m > 1e-10:
                c = np.where(abs(vd) > m * 1e-4)[0]
                if len(c) > 0:
                    if np.all(l[c] > 0):
                        vd[c] = vd[c] / l[c]
                        m = abs(vd[c]).max()
                    vd[c] /= m
                    if self.ont._is_s_begin:
                        plt = self.lattice.bar(
                            s[c] + l[c] / 2, vd[c], l[c], picker=True
                        )  # changed
                    else:
                        plt = self.lattice.bar(
                            s[c] - l[c] / 2, vd[c], l[c], picker=True
                        )  # changed
                    pl.setp(plt, facecolor=color, edgecolor=color)
                    if plt:
                        self.lines.append(plt[0])
                        self.legends.append(lbl)
                    row_names = self.ont.name
                    for r, i in zip(plt, c):
                        r.elemname = row_names[i]
                        r.elemprop = vdname
                        r.elemvalue = getattr(self.ont, vdname)[i]
                self.lattice.set_ylim(-1.5, 1.5)

    #    timeit('end lattice')

    def _column(self, name, sp, color):
        fig, s = self.figure, self.xaxis
        y = self.ont(name)[self.idx]
        (bxp,) = sp.plot(s, y, color, label=_mylbl(self.lglabel, name))
        sp.set_ylabel(_mylbl(self.axlabel, name))
        self.lines.append(bxp)
        self.legends.append(_mylbl(self.lglabel, name))
        sp.autoscale_view()

    def savefig(self, name):
        self.figure.savefig(name)
        return self


class Footprint(object):
    # class Footprint(ObjDebug):
    def __init__(self, x, y, tunx, tuny, nsigma, nangles, label="detuning"):
        self.nsigma = nsigma
        self.nangles = nangles
        self.x = np.array(x)
        self.y = np.array(y)
        self.tunx = np.array(tunx)
        self.tuny = np.array(tuny)
        self.label = label.replace("_", " ")

    def plot_grid(self, nsigma=None, lw=1):
        if nsigma is None:
            nsigma = self.nsigma
        nangles = self.nangles
        ranges = self.mkranges(nsigma)
        for i in ranges:
            if hasattr(i, "step"):
                lw = i.step == nangles and i.start / 2.0 or 1
            pl.plot(self.x[i], self.y[i], "-k", lw=lw)

    def mkranges(self, nsigma=None):
        if nsigma is None:
            nsigma = self.nsigma
        nangles = self.nangles
        ranges = []
        for i in range(nangles):
            ranges.append([0, i])
        for i in range(nangles):
            ranges.append(slice(1 + i, nangles * nsigma + 1, nangles))
        for i in range(nsigma):
            ranges.append(slice(1 + nangles * i, 1 + nangles * (i + 1)))
        return ranges

    def plot_footprint(
        t, nsigma=None, wp=(0.28, 0.31), spread=0.01, label=None, color=None
    ):
        ranges = t.mkranges(nsigma)
        lw = 1
        out = []
        lbl = True
        if label is None:
            label = t.label
        if color is None:
            color = "k"
        for i in ranges:
            if lbl:
                p = pl.plot(t.tunx[i], t.tuny[i], "-%s" % color, lw=lw, label=label)
                lbl = False
            else:
                p = pl.plot(t.tunx[i], t.tuny[i], "-%s" % color, lw=lw)
            out.append(p[0])
        pl.ylabel("$Q_y$")
        pl.xlabel("$Q_x$")
        pl.grid(True)
        qx, qy = wp
        pl.xlim(qx - spread, qx + spread)
        pl.ylim(qy - spread, qy + spread)
        return out

    def triangulate(t):
        tr = matplotlib.delaunay.triangulate.Triangulation(t.tunx, t.tuny)
        for i in tr.triangle_nodes:
            pl.plot(t.tunx[i], t.tuny[i])

    def reshape(self):
        """return tunes in [sigma,angles]"""
        qx = self.tunx[1:].reshape(self.nsigma, self.nangles)
        qy = self.tuny[1:].reshape(self.nsigma, self.nangles)
        return qx, qy


class FootTrack(Footprint):
    def __init__(self, dynapfn, nangles=7, nsigma=12, label="dynap"):
        self.label = label.replace("_", " ")
        t = tfsdata.open(dynapfn)
        self.tunx = t["tunx"]
        self.tuny = t["tuny"]
        self.tx = t["x"]
        self.ty = t["y"]
        self.nangles = nangles
        self.nsigma = nsigma
        # self.t=t


def mk_grid(nsigma, nangles):
    small = 0.05
    big = np.sqrt(1.0 - small**2)
    n = 1
    m = 0
    # sigma angle multiplier
    x = [small]
    y = [small]
    while n <= nsigma:
        angle = 90.0 / (nangles - 1) * m * np.pi / 180
        if m == 0:
            xs = n * big
            ys = n * small
        elif m == nangles - 1:
            xs = n * small
            ys = n * big
        else:
            xs = n * np.cos(angle)
            ys = n * np.sin(angle)
        m = m + 1
        if m == nangles:
            m = 0
            n = n + 1
        x.append(xs)
        y.append(ys)
    return np.array(x), np.array(y)


def Dq2(b2, Bx, Jx, By, Jy):
    #  b2=k1l
    Dqx = b2 * Bx / (4.0 * np.pi)
    Dqy = -b2 * By / (4.0 * np.pi)
    return Dqx, Dqy


def Dq4(b4, Bx, Jx, By, Jy):
    #  b4=k3l/6.
    Dqx = b4 * 3 * Bx * (Bx * Jx - 2 * By * Jy) / (8.0 * np.pi)
    Dqy = b4 * 3 * By * (By * Jy - 2 * Bx * Jx) / (8.0 * np.pi)
    return Dqx, Dqy


def Dq6(b6, Bx, Jx, By, Jy):
    #  b6=k5l/120.
    Dqx = (
        b6
        * 5
        * Bx
        * (Bx**2 * Jx**2 - 6 * Bx * By * Jx * Jy + 3 * By**2 * Jy**2)
        / (8.0 * np.pi)
    )
    Dqy = (
        -b6
        * 5
        * By
        * (3 * Bx**2 * Jx**2 - 6 * Bx * By * Jx * Jy + By**2 * Jy**2)
        / (8.0 * np.pi)
    )
    return Dqx, Dqy


def Dq8(b8, Bx, Jx, By, Jy):
    #  b8=k7l/5040.
    Dqx = (
        b8
        * 35
        * Bx
        * (
            Bx**3 * Jx**3
            - 12 * Bx**2 * By * Jx**2 * Jy
            + 18 * Bx * By**2 * Jx * Jy**2
            - 4 * By**3 * Jy**3
        )
        / (32.0 * np.pi)
    )
    Dqy = (
        b8
        * 35
        * By
        * (
            -4 * Bx**3 * Jx**3
            + 18 * Bx**2 * By * Jx**2 * Jy
            - 12 * Bx * By**2 * Jx * Jy**2
            + By**3 * Jy**3
        )
        / (32.0 * np.pi)
    )
    return Dqx, Dqy


def Dq10(b10, Bx, Jx, By, Jy):
    #  b10=k9l/362880.
    Dqx = (
        b10
        * 63
        * Bx
        * (
            Bx**4 * Jx**4
            - 20 * Bx**3 * By * Jx**3 * Jy
            + 60 * Bx**2 * By**2 * Jx**2 * Jy**2
            - 40 * Bx * By**3 * Jx * Jy**3
            + 5 * By**4 * Jy**4
        )
        / (32.0 * np.pi)
    )
    Dqy = (
        -b10
        * 63
        * By
        * (
            5 * Bx**4 * Jx**4
            - 40 * Bx**3 * By * Jx**3 * Jy
            + 60 * Bx**2 * By**2 * Jx**2 * Jy**2
            - 20 * Bx * By**3 * Jx * Jy**3
            + By**4 * Jy**4
        )
        / (32.0 * np.pi)
    )
    return Dqx, Dqy


def Dq12(b12, Bx, Jx, By, Jy):
    #  b12=k11l/39916800.
    Dqx = (
        b12
        * 231
        * Bx
        * (
            Bx**5 * Jx**5
            - 30 * Bx**4 * By * Jx**4 * Jy
            + 150 * Bx**3 * By**2 * Jx**3 * Jy**2
            - 200 * Bx**2 * By**3 * Jx**2 * Jy**3
            + 75 * Bx * By**4 * Jx * Jy**4
            - 6 * By**5 * Jy**5
        )
        / (64.0 * np.pi)
    )
    Dqy = (
        b12
        * 231
        * By
        * (
            -6 * Bx**5 * Jx**5
            + 75 * Bx**4 * By * Jx**4 * Jy
            - 200 * Bx**3 * By**2 * Jx**3 * Jy**2
            + 150 * Bx**2 * By**3 * Jx**2 * Jy**3
            - 30 * Bx * By**4 * Jx * Jy**4
            + By**5 * Jy**5
        )
        / (64.0 * np.pi)
    )
    return Dqx, Dqy


def Dq14(b14, Bx, Jx, By, Jy):
    #  b14=k13l/6227020800.
    Dqx = (
        b14
        * 429
        * Bx
        * (
            Bx**6 * Jx**6
            - 42 * Bx**5 * By * Jx**5 * Jy
            + 315 * Bx**4 * By**2 * Jx**4 * Jy**2
            - 700 * Bx**3 * By**3 * Jx**3 * Jy**3
            + 525 * Bx**2 * By**4 * Jx**2 * Jy**4
            - 126 * Bx * By**5 * Jx * Jy**5
            + 7 * By**6 * Jy**6
        )
        / (64.0 * np.pi)
    )
    Dqy = (
        -b14
        * 429
        * By
        * (
            7 * Bx**6 * Jx**6
            - 126 * Bx**5 * By * Jx**5 * Jy
            + 525 * Bx**4 * By**2 * Jx**4 * Jy**2
            - 700 * Bx**3 * By**3 * Jx**3 * Jy**3
            + 315 * Bx**2 * By**4 * Jx**2 * Jy**4
            - 42 * Bx * By**5 * Jx * Jy**5
            + By**6 * Jy**6
        )
        / (64.0 * np.pi)
    )
    return Dqx, Dqy


def Dq16(b16, Bx, Jx, By, Jy):
    #  b16=k15l/1307674368000.
    Dqx = (
        b16
        * 6435
        * Bx
        * (
            Bx**7 * Jx**7
            - 56 * Bx**6 * By * Jx**6 * Jy
            + 588 * Bx**5 * By**2 * Jx**5 * Jy**2
            - 1960 * Bx**4 * By**3 * Jx**4 * Jy**3
            + 2450 * Bx**3 * By**4 * Jx**3 * Jy**4
            - 1176 * Bx**2 * By**5 * Jx**2 * Jy**5
            + 196 * Bx * By**6 * Jx * Jy**6
            - 8 * By**7 * Jy**7
        )
        / (512.0 * np.pi)
    )
    Dqy = (
        b16
        * 6435
        * By
        * (
            -8 * Bx**7 * Jx**7
            + 196 * Bx**6 * By * Jx**6 * Jy
            - 1176 * Bx**5 * By**2 * Jx**5 * Jy**2
            + 2450 * Bx**4 * By**3 * Jx**4 * Jy**3
            - 1960 * Bx**3 * By**4 * Jx**3 * Jy**4
            + 588 * Bx**2 * By**5 * Jx**2 * Jy**5
            - 56 * Bx * By**6 * Jx * Jy**6
            + By**7 * Jy**7
        )
        / (512.0 * np.pi)
    )
    return Dqx, Dqy


def nchoosek(n, k):
    bc = [1 for i in range(0, k + 1)]
    for j in range(1, n - k + 1):
        for i in range(1, k + 1):
            bc[i] = bc[i - 1] + bc[i]
    return bc[k]


def factorial(n):
    fact = 1
    for x in range(1, n + 1):
        fact *= x
    return fact


def FeedDown(bn, an, x, y, i):
    cn = [b + 1j * a for b, a in zip(bn, an)]
    z = x + 1j * y
    n = len(cn)
    fd = sum([nchoosek(k, i) * cn[k - 1] * z ** (k - i) for k in range(i, n + 1)])
    return fd.real, fd.imag


def k2b(kn, ks, x=0, y=0):
    bn = [k / float(factorial(n)) for n, k in enumerate(kn)]
    an = [k / float(factorial(n)) for n, k in enumerate(ks)]
    z = complex(x, y)
    cn = [complex(b, a) for b, a in zip(bn, an)]
    n = len(cn)
    zn = [z**k for k in range(n)]
    fn = [
        sum([nchoosek(k, i) * cn[k - 1] * z ** (k - i) for k in range(i, n + 1)])
        for i in range(1, n + 1)
    ]
    return fn


def twiss2map(bet1, alf1, bet2, alf2, mu):
    b1b2 = np.sqrt(bet1 * bet2)
    b1onb2 = np.sqrt(bet1 / bet2)
    c = np.cos(2 * np.pi * mu)
    s = np.sin(2 * np.pi * mu)
    r11 = (c + alf1 * s) / b1onb2
    r12 = b1b2 * s
    r21 = ((alf1 - alf2) * c - (1 + alf1 * alf2) * s) / b1b2
    r22 = b1onb2 * (c - alf2 * s)
    return [[r11, r12], [r21, r22]]


mycolors = list("rcgmb")
