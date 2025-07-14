from __future__ import print_function

import time, os, re

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scipy.optimize

try:
    import pytimber
    import pjlsa
except ImportError:
    print("Warning: pytimber and/or pjlsa not installed")
from . import poly_fit

from .cubic_prep import cubic_fit1, cubic_fit2, cubic_val


def get_smooth(n, x, y):
    x0 = [x[0], x[-1]]
    y0 = [y[0], y[-1]]
    xp0 = [x[-1]]
    yp0 = [0]
    pol = poly_fit.poly_fit(n, x, y, x0, y0, xp0, yp0)
    return poly_fit.poly_val(pol, x)


def get_inequality(x, xmin, xmax, label=None, debug=False):
    xxmin = x.min()
    xxmax = x.max()
    pen = 0
    if xxmin < xmin:
        pen += xxmin**2
    if xxmax > xmax:
        pen += xxmax**2
    if debug:
        print("%s: %s<%s, %s>%s: %s" % (label, xxmin, xmin, xxmax, xmax, pen))
    return pen


clight = 299792458


def trim_range(vv, vmin, vmax):
    vv = vv.copy()
    vv[(vv > vmin) & (vv < vmax)] = 0
    return np.abs(vv)


def derivative(t, v):
    dt = np.diff(t)
    return t[:-1] + dt, np.diff(v) / dt


def madname_from_pcname(pc):
    name = ".".join(pc.split(".")[2:]).lower()
    if name.startswith("rcb"):
        return "a" + name[1:]
    else:
        return "k" + name[1:]


def get_calib_mad(lsa, madname):
    print(madname)
    pcname = lsa.findPCNameByMadStrength(madname)
    if pcname is not None:
        if pcname.split(".")[0] in ["RQD", "RQF", "RB"]:
            pcname = pcname + "B1"
        cal = lsa.getCalibration(pcname)
        if cal is not None:
            return [cal.field, cal.current]


def check_i_limit(i, imin, imax, label):
    if i > imax or i < imin:
        raise ValueError("%s=%g exceed limit of [%g,%g]" % (label, i, imin, imax))


class Circuit2in1(object):
    def __init__(self, name=None, **nargs):
        self.name = name
        self.name1 = name
        self.name2 = name.replace("B1", "B2")
        self.madname1 = madname_from_pcname(self.name1)
        self.madname2 = madname_from_pcname(self.name2)
        self.__dict__.update(nargs)

    def copy(self):
        other = self.__class__(self.name)
        for k, v in self.__dict__.items():
            other.__dict__[k] = v
        return other

    def smooth_old(self, step, stretch=1, kind="cubic"):
        pc = self.copy()
        newt = np.arange(pc.t[0], pc.t[-1], step)
        newi1 = scipy.interpolate.interp1d(pc.t, pc.i1, kind="cubic")(newt)
        newi2 = scipy.interpolate.interp1d(pc.t, pc.i2, kind="cubic")(newt)
        newsteps = {pc.madname1: newi1, pc.madname2: newi2}
        newt[-10:] = newt[-11] + step * np.arange(1, 11) ** stretch
        ((newt[-10:] - newt[-10]) ** stretch)
        pc.mk_simul(newt, newsteps)
        return pc

    def smooth(self, step, smooth=None):
        pc = self.copy()
        newt = np.arange(pc.t[0], pc.t[-1] + step, step)
        newi1 = scipy.interpolate.UnivariateSpline(pc.t, pc.i1, s=smooth)
        newi2 = scipy.interpolate.UnivariateSpline(pc.t, pc.i2, s=smooth)
        diff1 = abs(newi1(pc.t) - pc.i1).max()
        diff2 = abs(newi2(pc.t) - pc.i2).max()
        print("diff %s %s %s" % (pc.madname1, diff1, diff2))
        newsteps = {pc.madname1: newi1(newt), pc.madname2: newi2(newt)}
        pc.mk_simul(newt, newsteps)
        return pc

    def __repr__(self):
        tmp = """\
%s(%r,
%%s)""" % (
            self.__class__.__name__,
            self.name,
        )
        attrs = [
            "name1",
            "name2",
            "madname1",
            "madname2",
            "calib1",
            "calib2",
            "r1",
            "rc",
            "r2",
            "l1",
            "l2",
            "c1",
            "c2",
            "i1max",
            "i1min",
            "v1max",
            "v1min",
            "i1pmax",
            "i1pmin",
            "v1pmax",
            "v1pmin",
            "i1ppmax",
            "i1ppmin",
            "polarity1",
            "i2max",
            "i2min",
            "v2max",
            "v2min",
            "i2pmax",
            "i2pmin",
            "v2pmax",
            "v2pmin",
            "i2ppmax",
            "i2ppmin",
            "polarity2",
        ]
        out = []
        for att in attrs:
            if hasattr(self, att):
                out.append("    %s=%r" % (att, getattr(self, att)))
        return tmp % (",\n".join(out))

    def get_meas(self, t1, t2, db=None):
        if db is None:
            db = pytimber.LoggingDB()
        ldb = pytimber.PageStore("pytimber_cache.db", "pytimber_cache")
        v1n = self.name1 + ":V_MEAS"
        v2n = self.name2 + ":V_MEAS"
        i1n = self.name1 + ":I_MEAS"
        i2n = self.name2 + ":I_MEAS"
        data = db.get([v1n, v2n, i1n, i2n], t1, t2)
        ldb.store(data)
        tv1, v1 = data[v1n]
        tv2, v2 = data[v2n]
        ti1, i1 = data[i1n]
        ti2, i2 = data[i2n]
        self.t = np.arange(tv1[0], tv1[-1], 0.2)
        self.v1 = np.interp(self.t, tv1, v1)
        self.v2 = np.interp(self.t, tv2, v2)
        self.i1 = np.interp(self.t, ti1, i1)
        self.i2 = np.interp(self.t, ti2, i2)
        self.i1p = np.interp(self.t, *derivative(ti1, i1))
        self.i2p = np.interp(self.t, *derivative(ti2, i2))
        self.v1p = np.interp(self.t, *derivative(tv1, v1))
        self.v2p = np.interp(self.t, *derivative(tv2, v2))
        self.i1pp = np.interp(self.t, *derivative(self.t, self.i1p))
        self.i2pp = np.interp(self.t, *derivative(self.t, self.i2p))
        return self

    def mk_simul(self, tsteps, isteps):
        i1 = isteps[self.madname1]
        i2 = isteps[self.madname2]
        self.t = np.array(tsteps)
        self.i1 = np.array(i1)
        self.i2 = np.array(i2)
        self.i1p = np.interp(self.t, *derivative(self.t, i1))
        self.i2p = np.interp(self.t, *derivative(self.t, i2))
        self.v1, self.v2 = self.model_getv()
        self.v1p = np.interp(self.t, *derivative(self.t, self.v1))
        self.v2p = np.interp(self.t, *derivative(self.t, self.v2))
        self.i1pp = np.interp(self.t, *derivative(self.t, self.i1p))
        self.i2pp = np.interp(self.t, *derivative(self.t, self.i2p))
        return self

    def mk_simul_summ_vec(self):
        ti1 = trim_range(self.i1, self.i1min, self.i1max)
        ti2 = trim_range(self.i2, self.i2min, self.i2max)
        ti1p = trim_range(self.i1p, self.i1pmin, self.i1pmax)
        ti2p = trim_range(self.i2p, self.i2pmin, self.i2pmax)
        ti1pp = trim_range(self.i1pp, self.i1ppmin, self.i1ppmax)
        ti2pp = trim_range(self.i2pp, self.i2ppmin, self.i2ppmax)
        tv1 = trim_range(self.v1, self.v1min, self.v1max)
        tv2 = trim_range(self.v2, self.v2min, self.v2max)
        tv1p = trim_range(self.v1p, self.v1pmin, self.v1pmax)
        tv2p = trim_range(self.v2p, self.v2pmin, self.v2pmax)
        summ1 = ti1 + ti1 + ti1p + tv1 + tv1p + ti1pp
        summ2 = ti2 + ti2 + ti2p + tv2 + tv2p + ti2pp
        return summ1 + summ2

    def mk_simul_summ(self):
        i1min = self.i1.min()
        i2min = self.i2.min()
        i1max = self.i1.max()
        i2max = self.i2.max()
        i1pmax = abs(self.i1p).max()
        i2pmax = abs(self.i2p).max()
        v1min = self.v1.min()
        v2min = self.v2.min()
        v1pmax = abs(self.v1p).max()
        v2pmax = abs(self.v2p).max()
        i1ppmax = abs(self.i1pp).max()
        i2ppmax = abs(self.i2pp).max()
        di1min = i1min - self.i1min
        di2min = i2min - self.i2min
        di1max = self.i1max - i1max
        di2max = self.i2max - i2max
        di1pmax = self.i1pmax - i1pmax
        di2pmax = self.i2pmax - i2pmax
        dv1min = v1min - self.v1min
        dv2min = v2min - self.v2min
        dv1pmax = self.v1pmax - v1pmax
        dv2pmax = self.v2pmax - v2pmax
        di1ppmax = self.i1ppmax - i1ppmax
        di2ppmax = self.i2ppmax - i2ppmax
        ti1min = 0 if di1min > 0 else di1min
        ti2min = 0 if di2min > 0 else di2min
        ti1max = 0 if di1max > 0 else di1max
        ti2max = 0 if di2max > 0 else di2max
        ti1pmax = 0 if di1pmax > 0 else di1pmax
        ti2pmax = 0 if di2pmax > 0 else di2pmax
        tv1min = 0 if dv1min > 0 else dv1min
        tv2min = 0 if dv2min > 0 else dv2min
        tv1pmax = 0 if dv1pmax > 0 else dv1pmax
        tv2pmax = 0 if dv2pmax > 0 else dv2pmax
        ti1ppmax = 0 if di1ppmax > 0 else di1ppmax
        ti2ppmax = 0 if di2ppmax > 0 else di2ppmax
        summ1 = ti1min + ti1max + ti1pmax + tv1min + tv1pmax + ti1ppmax
        summ2 = ti2min + ti2max + ti2pmax + tv2min + tv2pmax + ti2ppmax
        fmt = "%-15s %7.1f %7.1f %7.3f %7.3f %7.3f %7.3f %7.3f"
        if abs(summ1) > 0:
            print(
                fmt
                % (self.madname1, i1min, i1max, i1pmax, i1ppmax, v1min, v1pmax, summ1)
            )
        if abs(summ2) > 0:
            print(
                fmt
                % (self.madname2, i2min, i2max, i2pmax, i2ppmax, v2min, v2pmax, summ2)
            )
        return summ1 + summ2

    def plot_vi(self, num=None, vp=False, label1=None, label2=None):
        if num is None:
            num = self.name1
        if label1 is None:
            label1 = self.name1
        if label2 is None:
            label2 = self.name2
        f = pl.figure(num=num, figsize=(8, 8))
        f.clf()
        if vp:
            f, (ax1, ax2, ax3, ax4, ax5) = pl.subplots(5, sharex=True, num=num)
        else:
            f, (ax1, ax2, ax3, ax4) = pl.subplots(4, sharex=True, num=num)
        print(num)
        ax1.plot(self.t, self.i1, label=label1)
        ax1.plot(self.t, self.i2, label=label2)
        # ax1.set_xlabel('Time [s]')
        ax1.set_ylabel("Current [A]")
        ax1.legend()
        ax1.axhline(self.i1min, color="k")
        ax1.axhline(self.i1max, color="k")
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()
        ax2.plot(self.t, self.i1p, label=label1)
        ax2.plot(self.t, self.i2p, label=label2)
        # ax2.set_xlabel('Time [s]')
        ax2.set_ylabel("$I'(t)$  [A/s]")
        ax2.axhline(self.i1pmin, color="k")
        ax2.axhline(self.i1pmax, color="k")
        # ax2.legend()
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()
        ax3.plot(self.t, self.i1pp, label=label1)
        ax3.plot(self.t, self.i2pp, label=label2)
        # ax3.set_xlabel('Time [s]')
        ax3.set_ylabel(r"I''(t) [I/$\rm s^2$]")
        ax3.axhline(self.i1ppmin, color="k")
        ax3.axhline(self.i1ppmax, color="k")
        # ax3.legend(
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()
        ax4.plot(self.t, self.v1, label=label1)
        ax4.plot(self.t, self.v2, label=label2)
        ax4.set_xlabel("Time [s]")
        ax4.set_ylabel("$V(t)$ [V]")
        if self.v1min == 0:
            ax4.axhline(self.v1min, color="k")
        # ax4.axhline(self.v1max,color='k')
        # ax4.legend()
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()
        if vp:
            ax5.plot(self.t, self.v1p * 1e3, label=label1)
            ax5.plot(self.t, self.v2p * 1e3, label=label2)
            # ax5.set_xlabel('Time [s]')
            ax5.set_ylabel("$V'(t)$ [mV/s]")
            ax5.axhline(self.v1pmin * 1e3, color="k")
            ax5.axhline(self.v1pmax * 1e3, color="k")
            # ax5.legend()
        return self

    def model_set(self, x):
        self.r1, self.rc, self.r2, self.l1, self.l2, self.c1, self.c2 = x
        return self

    def model_get(self):
        return self.r1, self.rc, self.r2, self.l1, self.l2, self.c1, self.c2

    def model_getv(self):
        r1, rc, r2, l1, l2 = self.r1, self.rc, self.r2, self.l1, self.l2
        c1, c2 = self.c1, self.c2
        i1 = self.i1
        i2 = self.i2
        i1p = self.i1p
        i2p = self.i2p
        di = i1 - i2
        v1 = r1 * i1 + rc * di + l1 * i1p + c1
        v2 = r2 * i2 - rc * di + l2 * i2p + c2
        return v1, v2

    def model_get_constraints_cubic(
        self, cubic1, cubic2, tstep, sup=np.linspace(0, 1.0, 3), debug=False
    ):
        # specs
        if debug:
            print("constraints %s, tstep %s" % (self.madname1, tstep))
        i1min, i1pmin, i1ppmin = self.i1min, self.i1pmin, self.i1ppmin
        i1max, i1pmax, i1ppmax = self.i1max, self.i1pmax, self.i1ppmax
        v1min, v1pmin = self.v1min + self.c1, self.v1pmin
        v1max, v1pmax = self.v1max, self.v1pmax
        i2min, i2pmin, i2ppmin = self.i2min, self.i2pmin, self.i2ppmin
        i2max, i2pmax, i2ppmax = self.i2max, self.i2pmax, self.i2ppmax
        v2min, v2pmin = self.v2min + self.c2, self.v2pmin
        v2max, v2pmax = self.v2max, self.v2pmax
        r1, rc, r2, l1, l2 = self.r1, self.rc, self.r2, self.l1, self.l2
        c1, c2 = self.c1, self.c2
        # model
        t = tstep * sup
        i1, i1p, i1pp = cubic_val(t, *cubic1)
        i2, i2p, i2pp = cubic_val(t, *cubic2)
        di = i1 - i2
        v1 = r1 * i1 + rc * di + l1 * i1p + c1
        v2 = r2 * i2 - rc * di + l2 * i2p + c2
        pen = get_inequality(i1, i1min, i1max, "i1", debug=debug)
        pen += get_inequality(v1, v1min, v1max, "v1", debug=debug)
        pen += get_inequality(i1p, i1pmin, i1pmax, "i1p", debug=debug)
        pen += get_inequality(i1pp, i1ppmin, i1ppmax, "i1pp", debug=debug)
        pen += get_inequality(i2, i2min, i2max, "i2", debug=debug)
        pen += get_inequality(v2, v2min, v2max, "v2", debug=debug)
        pen += get_inequality(i2p, i2pmin, i2pmax, "i2p", debug=debug)
        pen += get_inequality(i2pp, i2ppmin, i2ppmax, "i2pp", debug=debug)
        return pen

    def model_res(self, x):
        self.model_set(x)
        v1mod, v2mod = self.model_getv()
        res1 = self.v1 - v1mod
        res2 = self.v2 - v2mod
        return res1**2 + res2**2

    def model_fit(self, x0=[1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0, 0]):
        res = scipy.optimize.leastsq(self.model_res, x0, full_output=True)
        x, cov, info, mesg, ier = res
        # print mesg
        self.model_set(x)
        return self

    def model_plot(self):
        v1mod, v2mod = self.model_getv()
        pl.plot(self.t, self.v1, label="V1 meas.")
        pl.plot(self.t, self.v2, label="V2 meas.")
        pl.plot(self.t, v1mod, label="V1 model")
        pl.plot(self.t, v2mod, label="V2 model")
        pl.xlabel("Time [s]")
        pl.ylabel("Voltage [V]")
        pl.legend()
        if self.t[0] > 1e6:
            pytimber.set_xaxis_date()

    def model_plot_res(self, xaxis_date=True):
        v1mod, v2mod = self.model_getv()
        pl.plot(self.t, self.v1, label=r"$V_1$")
        pl.plot(self.t, self.v2, label=r"$V_2$")
        pl.plot(self.t, self.v1 - v1mod, label=r"$\delta V_1$")
        pl.plot(self.t, self.v2 - v2mod, label=r"$\delta V_2$")
        pl.xlabel("Time [s]")
        pl.ylabel("Voltage [V]")
        pl.legend()
        if self.t[0] > 1e6 and xaxis_date is True:
            pytimber.set_xaxis_date()
        return self

    def set_calib_from_lsa(self, lsa=None):
        if lsa is None:
            lsa = pjlsa.LSAClient()
        self.calib1 = get_calib_mad(lsa, self.madname1)
        self.calib2 = get_calib_mad(lsa, self.madname2)

    def k2i(self, k1, k2, momentum):
        field1 = self.polarity1 * k1 * momentum / 299792458
        field2 = self.polarity2 * k2 * momentum / 299792458
        field1c, i1c = self.calib1
        field2c, i2c = self.calib2
        i1 = np.interp(field1, field1c, i1c, right=i1c[-1] / field1c[-1] * field1)
        i2 = np.interp(field2, field2c, i2c, right=i2c[-1] / field2c[-1] * field2)
        return i1, i2

    def i2field(self, i1, i2):
        field1c, i1c = self.calib1
        field2c, i2c = self.calib2
        field1 = np.interp(i1, i1c, field1c, right=i1 / i1c[-1] * field1c[-1])
        field2 = np.interp(i2, i2c, field2c, right=i2 / i2c[-1] * field2c[-1])
        return field1, field2

    def get_min_tstep(
        self,
        momentum,
        tstart,
        k1start,
        k1next,
        i1pstart,
        k2start,
        k2next,
        i2pstart,
        i1pforce=None,
        i2pforce=None,
        startstep=0,
        tstepmin=0.02,
    ):
        # specs
        i1min, i1pmin, i1ppmin = self.i1min, self.i1pmin, self.i1ppmin
        i1max, i1pmax, i1ppmax = self.i1max, self.i1pmax, self.i1ppmax
        v1min, v1pmin = self.v1min + self.c1, self.v1pmin
        v1max, v1pmax = self.v1max, self.v1pmax
        i2min, i2pmin, i2ppmin = self.i2min, self.i2pmin, self.i2ppmin
        i2max, i2pmax, i2ppmax = self.i2max, self.i2pmax, self.i2ppmax
        v2min, v2pmin = self.v2min + self.c2, self.v2pmin
        v2max, v2pmax = self.v2max, self.v2pmax
        if i1pforce is not None:
            i1pmin = -i1pforce
            i1pmax = i1pforce
        if i2pforce is not None:
            i2pmin = -i2pforce
            i2pmax = i2pforce
        r1, l1, c1 = self.r1, self.l1, self.c1
        r2, l2, c2 = self.r2, self.l2, self.c2
        # endspecs
        # initial conditions
        i1start, i2start = self.k2i(k1start, k2start, momentum(tstart))
        check_i_limit(i1start, i1min, i1max, self.name)
        check_i_limit(i2start, i2min, i2max, self.name)
        di = i1start - i2start
        rc = self.rc
        v1start = c1 + r1 * i1start + l1 * i1pstart + rc * di
        v2start = c2 + r2 * i2start + l2 * i2pstart - rc * di

        def ftomin(x):
            (tstep,) = x
            i1next, i2next = self.k2i(k1next, k2next, momentum(tstart + tstep))
            i1p = (i1next - i1start) / tstep
            i2p = (i2next - i2start) / tstep
            i1pp = (i1p - i1pstart) / tstep
            i2pp = (i2p - i2pstart) / tstep
            v1next = c1 + r1 * i1next + l1 * i1p + rc * di
            v2next = c2 + r2 * i2next + l2 * i2p - rc * di
            v1p = (v1next - v1start) / tstep
            v2p = (v2next - v2start) / tstep
            # pen=tstep**2*1e-9
            # if tstep < 0.1:
            #   print(tstep)
            pen = 0
            if i1next > i1max or i1next < i1min:
                pen += i1next**2
            if i2next > i2max or i2next < i2min:
                pen += i2next**2
            if i1p > i1pmax or i1p < i1pmin:
                pen += i1p**2
            #              print('i1p',i1p)
            if i2p > i2pmax or i2p < i2pmin:
                pen += i2p**2
            #              print('i2p',i2p)
            if i1pp > i1ppmax or i1pp < i1ppmin:
                pen += i1pp**2
            #              print('i1p',i1p)
            if i2pp > i2ppmax or i2pp < i2ppmin:
                pen += i2pp**2
            #              print('i2p',i2p)
            if v1p > v1pmax or v1p < v1pmin:
                pen += v1p**2
            #      print('v1p',v1p)
            if v2p > v2pmax or v2p < v2pmin:
                pen += v2p**2
            #    if tstep < 0.1:
            #      print('v2p',v2p) if tstep < 0.1
            if v1next > v1max or v1next < v1min:
                pen += v1next**2
            #             print('v1next',v1next)
            if v2next > v2max or v2next < v2min:
                pen += v2next**2
            #              print('v2next',v2next)
            # if pen==0:
            #  pen=tstep**2*1e-9
            # print('pen',pen)
            # if (tstart>570 and tstart < 590 ):
            #  print(tstart,tstep,v1start,v1next,v1p)
            return pen

        # tstep,=scipy.optimize.fmin(ftomin,1e-9,ftol=1e-16,
        #          xtol=1e-9,disp=False)
        # linear scan
        tstep = startstep
        pen = ftomin([tstep])
        while tstep < 20.0 and pen > 0:
            tstep += tstepmin
            pen = ftomin([tstep])
        ## bisection scan
        # tstep=startstep
        # pen=ftomin([tstep])
        # if pen>0:
        #  t1=startstep
        #  t2=1000
        #  while (t2-t1)>1e-6:
        #    t3=(t1+t2)/2
        #    pen=ftomin([t3])
        #    if pen==0:
        #        t2=t3
        #    else:
        #        t1=t3
        #  tstep=t3*1.2
        ## end bisection scan
        i1next, i2next = self.k2i(k1next, k2next, momentum(tstart + tstep))
        check_i_limit(i1next, i1min, i1max, self.name)
        check_i_limit(i2next, i2min, i2max, self.name)
        if tstep <= 1e-9:
            return 0, 0, 0
        else:
            i1pnext = (i1next - i1start) / tstep
            i2pnext = (i2next - i2start) / tstep
            return tstep, i1pnext, i2pnext

    def get_min_tstep2(
        self, momentum, tstart, k1steps, k2steps, startstep, tstepmin=0.02, maxstep=20
    ):
        # specs
        i1min, i1pmin, i1ppmin = self.i1min, self.i1pmin, self.i1ppmin
        i1max, i1pmax, i1ppmax = self.i1max, self.i1pmax, self.i1ppmax
        v1min, v1pmin = self.v1min + self.c1, self.v1pmin
        v1max, v1pmax = self.v1max, self.v1pmax
        i2min, i2pmin, i2ppmin = self.i2min, self.i2pmin, self.i2ppmin
        i2max, i2pmax, i2ppmax = self.i2max, self.i2pmax, self.i2ppmax
        v2min, v2pmin = self.v2min + self.c2, self.v2pmin
        v2max, v2pmax = self.v2max, self.v2pmax
        r1, l1, c1 = self.r1, self.l1, self.c1
        r2, l2, c2 = self.r2, self.l2, self.c2
        rc = self.rc
        # endspecs
        # initial conditions
        t_2 = tstart
        k1_0 = k1steps[-1]
        k1_1 = k1steps[-2]
        k1_2 = k1steps[-3]
        k2_0 = k2steps[-1]
        k2_1 = k2steps[-2]
        k2_2 = k2steps[-3]

        # i1start,i2start=self.k2i(k1start, k2start, momentum(tstart))
        # check_i_limit(i1start,i1min,i1max,self.name)
        # check_i_limit(i2start,i2min,i2max,self.name)
        # di=i1start-i2start
        # v1start=c1+r1*i1start+l1*i1pstart+rc*di
        # v2start=c2+r2*i2start+l2*i2pstart-rc*di
        def ftomin(x):
            (tstep,) = x
            t_1 = t_2 + tstep
            t_0 = t_1 + tstep
            i1_0, i2_0 = self.k2i(k1_0, k2_0, momentum(t_0))
            i1_1, i2_1 = self.k2i(k1_1, k2_1, momentum(t_1))
            i1_2, i2_2 = self.k2i(k1_2, k2_2, momentum(t_2))
            i1p = (i1_0 - i1_1) / tstep
            i2p = (i2_0 - i2_1) / tstep
            i1pp = (i1_0 - 2 * i1_1 + i1_2) / tstep**2
            i2pp = (i2_0 - 2 * i2_1 + i2_2) / tstep**2
            di = i1_1 - i2_1
            dip = i1p - i2p
            v1 = c1 + r1 * i1_1 + l1 * i1p + rc * di
            v2 = c2 + r2 * i2_1 + l2 * i2p - rc * di
            v1p = c1 + r1 * i1p + l1 * i1pp + rc * dip
            v2p = c2 + r2 * i2p + l2 * i2pp - rc * dip
            # pen=tstep**2*1e-9
            # if tstep < 0.1:
            #   print(tstep)
            pen = 0
            if i1_0 > i1max or i1_0 < i1min:
                pen += i1_0**2
                # print("pen %-30s: i1  %g"%(self.madname1,i1_0))
            if i2_0 > i2max or i2_0 < i2min:
                pen += i2_0**2
                # print("pen %-30s: i2  %g"%(self.madname2,i2_0))
            if i1p > i1pmax or i1p < i1pmin:
                pen += i1p**2
                # print("pen %-30s: i1p  %g"%(self.madname1,i1p))
            if i2p > i2pmax or i2p < i2pmin:
                pen += i2p**2
                # print("pen %-30s: i2p  %g"%(self.madname2,i2p))
            if i1pp > i1ppmax or i1pp < i1ppmin:
                pen += i1pp**2
                # print("pen %-30s: i1pp  %g"%(self.madname1,i1pp))
            if i2pp > i2ppmax or i2pp < i2ppmin:
                pen += i2pp**2
                # print("pen %-30s: i2pp  %g"%(self.madname2,i2pp))
            if v1p > v1pmax or v1p < v1pmin:
                pen += v1p**2
                # print("pen %-30s: v1p  %g"%(self.madname1,v1p))
            if v2p > v2pmax or v2p < v2pmin:
                pen += v2p**2
                # print("pen %-30s: v2p  %g"%(self.madname1,v2p))
            #    if tstep < 0.1:
            #      print('v2p',v2p) if tstep < 0.1
            if v1 > v1max or v1 < v1min:
                pen += v1**2
                # print("pen %-30s: v1  %g"%(self.madname1,v1))
            #             print('v1next',v1next)
            if v2 > v2max or v2 < v2min:
                pen += v2**2
                # print("pen %-30s: v2  %g"%(self.madname2,v2))
            #              print('v2next',v2next)
            # if pen==0:
            #  pen=tstep**2*1e-9
            # print('pen',pen)
            # if (tstart>570 and tstart < 590 ):
            #  print(tstart,tstep,v1start,v1next,v1p)
            return pen

        # tstep,=scipy.optimize.fmin(ftomin,1e-9,ftol=1e-16,
        #          xtol=1e-9,disp=False)
        # linear scan
        tstep = startstep
        pen = ftomin([tstep])
        while tstep < maxstep and pen > 0:
            tstep += tstepmin
            pen = ftomin([tstep])
        ## bisection scan
        # tstep=startstep
        # pen=ftomin([tstep])
        # if pen>0:
        #  t1=startstep
        #  t2=1000
        #  while (t2-t1)>1e-6:
        #    t3=(t1+t2)/2
        #    pen=ftomin([t3])
        #    if pen==0:
        #        t2=t3
        #    else:
        #        t1=t3
        #  tstep=t3*1.2
        ## end bisection scan
        # check_i_limit(i1next,i1min,i1max,self.name)
        # check_i_limit(i2next,i2min,i2max,self.name)
        if tstep <= 1e-9:
            return 0
        else:
            return tstep

    def get_min_kstep(
        self,
        pcstart,
        pcnext,
        k1func,
        k2func,
        kstart,
        kstop,
        tstep,
        i1pstart,
        i2pstart,
        i1pforce=None,
        i2pforce=None,
    ):
        # specs
        i1min, i1pmin = self.i1min, self.i1pmin
        i1max, i1pmax = self.i1max, self.i1pmax
        v1min, v1pmin = self.v1min + self.c1, self.v1pmin
        v1max, v1pmax = self.v1max, self.v1pmax
        i2min, i2pmin = self.i2min, self.i2pmin
        i2max, i2pmax = self.i2max, self.i2pmax
        v2min, v2pmin = self.v2min + self.c2, self.v2pmin
        v2max, v2pmax = self.v2max, self.v2pmax
        if i1pforce is not None:
            i1pmin = -i1pforce
            i1pmax = i1pforce
        if i2pforce is not None:
            i2pmin = -i2pforce
            i2pmax = i2pforce
        r1, l1, c1 = self.r1, self.l1, self.c1
        r2, l2, c2 = self.r2, self.l2, self.c2
        # endspecs
        # initial conditions
        k1start = k1func(kstart)
        k2start = k2func(kstart)
        i1start, i2start = self.k2i(k1start, k2start, pcstart)
        check_i_limit(i1start, i1min, i1max, self.name)
        check_i_limit(i2start, i2min, i2max, self.name)
        di = i1start - i2start
        rc = self.rc
        v1start = c1 + r1 * i1start + l1 * i1pstart + rc * di
        v2start = c2 + r2 * i2start + l2 * i2pstart - rc * di

        def ftomin(x):
            (knext,) = x
            k1next = k1func(knext)
            k2next = k2func(knext)
            i1next, i2next = self.k2i(k1next, k2next, pcnext)
            i1p = (i1next - i1start) / tstep
            i2p = (i2next - i2start) / tstep
            v1next = c1 + r1 * i1next + l1 * i1p + rc * di
            v2next = c2 + r2 * i2next + l2 * i2p - rc * di
            v1p = (v1next - v1start) / tstep
            v2p = (v2next - v2start) / tstep
            pen = 0
            if i1next > i1max or i1next < i1min:
                pen += i1next**2
                print("pen %-30s: i1  %g" % (self.madname1, i1next))
            if i1p > i1pmax or i1p < i1pmin:
                pen += i1p**2
                print("pen %-30s: i1p %g" % (self.madname1, i1p))
            if v1next > v1max or v1next < v1min:
                pen += v1next**2
                print("pen %-30s: v1  %g" % (self.madname1, v1next))
            if v1p > v1pmax or v1p < v1pmin:
                pen += v1p**2
                print("pen %-30s: v1p %g" % (self.madname1, v1p))
            if i2next > i2max or i2next < i2min:
                pen += i2next**2
                print("pen %-30s: i2  %g" % (self.madname2, i2next))
            if i2p > i2pmax or i2p < i2pmin:
                pen += i2p**2
                print("pen %-30s: ip2 %g" % (self.madname2, i2p))
            if v2next > v2max or v2next < v2min:
                pen += v2next**2
                print("pen %-30s: v2  %g" % (self.madname2, v2next))
            if v2p > v2pmax or v2p < v2pmin:
                pen += v2p**2
                # print("pen %-30s: vp2 %g"%(self.madname2,v2p))
            return pen

        knext = scipy.optimize.fmin_l_bfgs_b(
            ftomin, [kstart], approx_grad=True, disp=False, bounds=[(kstart, kstop)]
        )
        pen = ftomin([knext])
        if pen > 0:
            raise ValueError
        dist = abs(kstop - knext)
        kstep = dist / 16.0
        eps = dist * 1e-10
        while knext < kstop and pen == 0:
            knext += kstep
            pen = ftomin([knext])
            nnn += 1
        print(nnn, knext, pen, kstop)
        if kstop - knext < eps and pen == 0:  # always good
            knext = kstop
        else:  # now branching
            k1 = knext - kstep
            k2 = knext
            print("Branching", k1, k2, kstop)
            while abs(k2 - k1) > eps and nnn < 100:
                k3 = (k2 + k1) / 2
                pen = ftomin([k3])
                # print(("%15.7g"*5)%(pen,k1,k3,k2,k2-k1))
                if pen == 0:
                    k1 = k3
                else:
                    k2 = k3
                nnn += 1
            knext = k1
            if knext >= kstop:
                knext = kstop
        pen = ftomin([knext])
        if pen > 0:  # reduce step size
            print("Penalty", self.madname1, knext, pen)
            t1 = 0
            t2 = tstep
        fmt = "Opt %-20s %12.5g %12.5g %12.5g %4d"
        print(fmt % (self.madname1, kstart, kstop, knext - kstart, nnn))
        k1next = k1func(knext)
        k2next = k2func(knext)
        i1next, i2next = self.k2i(k1next, k2next, pcnext)
        check_i_limit(i1next, i1min, i1max, self.name)
        check_i_limit(i2next, i2min, i2max, self.name)
        i1pnext = (i1next - i1start) / tstep
        i2pnext = (i2next - i2start) / tstep
        return knext, i1pnext, i2pnext

    def set_pcinfo_from_lsa(self, lsa=None):
        if lsa is None:
            lsa = pjlsa.LSAClient()
        pc = lsa.getPCInfo(self.name1)
        self.i1pmin = pc.didtMin
        self.i1pmax = pc.didtMax
        self.i1ppmin = pc.decelerationLimit
        self.i1ppmax = pc.accelerationLimit
        self.i1min = pc.iMinOp
        self.i1max = pc.iNom
        # self.i1max=pc.iPNo
        # self.i1nom=pc.iNom
        # self.i1ult=pc.iUlt
        # self.polarity1=pc.polaritySwitch
        pc = lsa.getPCInfo(self.name2)
        self.i2pmin = pc.didtMin
        self.i2pmax = pc.didtMax
        self.i2min = pc.iMinOp
        self.i2max = pc.iNom
        self.i2ppmin = pc.decelerationLimit
        self.i2ppmax = pc.accelerationLimit
        # self.i2max=pc.iPNo
        # self.i2nom=pc.iNom
        # self.i2ult=pc.iUlt
        # self.polarity2=pc.polaritySwitch

    def set_vprop(self, vmin, vmax, vpmin, vpmax):
        self.v1max = vmax
        self.v1pmax = vpmax
        self.v1min = vmin
        self.v1pmin = vpmin
        self.v2max = vmax
        self.v2pmax = vpmax
        self.v2min = vmin
        self.v2pmin = vpmin
        if vmin < 0:
            self.i1min = -self.i1max
            self.i2min = -self.i2max


class CircuitSingle(object):
    def __init__(self, name, **nargs):
        self.name = name
        self.madname = madname_from_pcname(name)
        self.__dict__.update(nargs)

    def copy(self):
        other = self.__class__(self.name)
        for k, v in self.__dict__.items():
            other.__dict__[k] = v
        return other

    def __repr__(self):
        tmp = """\
%s(%r,
%%s)""" % (
            self.__class__.__name__,
            self.name,
        )
        attrs = [
            "madname",
            "calib",
            "r1",
            "l1",
            "c1",
            "imax",
            "imin",
            "vmax",
            "vmin",
            "ipmax",
            "ipmin",
            "vpmax",
            "vpmin",
            "polarity1",
        ]
        ldb = pytimber.PageStore("pytimber_cache.db", "pytimber_cache")
        out = []
        for att in attrs:
            if hasattr(self, att):
                out.append("     %s=%r" % (att, getattr(self, att)))
        return tmp % (",\n".join(out))

    def get_meas(self, t1, t2, db=None):
        if db is None:
            db = pytimber.LoggingDB()
        ldb = pytimber.PageStore("pytimber_cache.db", "pytimber_cache")
        v1n = self.name + ":V_MEAS"
        i1n = self.name + ":I_MEAS"
        data = db.get([v1n, i1n], t1, t2)
        ldb.store(data)
        tv1, v1 = data[v1n]
        ti1, i1 = data[i1n]
        self.t = np.arange(tv1[0], tv1[-1], 0.2)
        self.v1 = np.interp(self.t, tv1, v1)
        self.i1 = np.interp(self.t, ti1, i1)
        self.i1p = np.interp(self.t, *derivative(ti1, i1))
        self.v1p = np.interp(self.t, *derivative(tv1, v1))
        self.i1pp = np.interp(self.t, *derivative(self.t, self.i1p))

        return self

    def k2i(self, k1, momentum):
        field1 = self.polarity1 * k1 * momentum / 299792458
        field1c, i1c = self.calib
        i1 = np.interp(field1, field1c, i1c, right=field1 * i1c[-1] / field1c[-1])
        return i1

    def plot_vi(self):
        pl.subplot(211)
        pl.plot(self.t, self.v1)
        pytimber.set_xaxis_date()
        pl.subplot(212)
        pl.plot(self.t, self.i1)
        pytimber.set_xaxis_date()
        return self

    def model_set(self, x):
        self.r1, self.l1, self.c1 = x
        return self

    def model_get(self):
        return self.r1, self.l1, self.c1

    def model_getv(self):
        r1, l1, c1 = self.r1, self.l1, self.c1
        i1 = self.i1
        i1p = self.i1p
        v1 = r1 * i1 + l1 * i1p + c1
        return v1

    def model_res(self, x):
        self.model_set(x)
        v1mod = self.model_getv()
        res1 = self.v1 - v1mod
        return res1

    def model_fit(self, x0=[1e-3, 1e-3, 0]):
        res = scipy.optimize.leastsq(self.model_res, x0, full_output=True)
        x, cov, info, mesg, ier = res
        # self.model_fit_info=(info,mesg,ier)
        # print mesg
        self.model_set(x)
        return self

    def model_plot(self):
        v1mod = self.model_getv()
        pl.plot(self.t, self.v1, label="V1 meas.")
        pl.plot(self.t, v1mod, label="V1 model")
        pytimber.set_xaxis_date()
        return self

    def set_calib_from_lsa(self, lsa=None):
        if lsa is None:
            lsa = pjlsa.LSAClient()
        self.calib1 = get_calib_mad(lsa, self.madname)

    def set_pcinfo_from_lsa(self, lsa=None):
        if lsa is None:
            lsa = pjlsa.LSAClient()
        pc = lsa.getPCInfo(self.name)
        self.ipmin = pc.didtMin
        self.ipmax = pc.didtMax
        self.imin = pc.iMinOp
        self.imax = pc.iNom
        # self.imax=pc.iPNo
        # self.inom=pc.iNom
        # self.iult=pc.iUlt
        # self.polarity=pc.polaritySwitch

    def get_min_tstep(self, momentum, tstart, k1start, k1next, i1pstart, i1start=None):
        i1next = self.k2i(k1next, momentum(tnext))
        if i1next > i1max or i1next < i1min:
            raise ValueError("I=%g exceed limit of [%g,%g]" % (i1next, i1min, i1max))
        if i1start is not None:
            i1start = self.k2i(k1start, momentum(tstart))
        v1start = c1 + r1 * i1start + l1 * i1pstart
        r1, l1, c1 = self.r1, self.l1, self.c1
        i1max, i1pmax = self.i1max, self.i1pmax
        v1max, v1pmax = self.v1max, self.v1pmax
        i1min, i1pmin = self.i1min, self.i1pmin
        v1min, v1pmin = self.v1min, self.v1pmin

        def tomin(tstep):
            i1next = self.k2i(k1next, momentum(tstart + tstep))
            i1p = (i1next - i1start) / tstep
            v1next = c1 + r1 * i1next + l1 * i1p
            v1p = (v1next - v1start) / tstep
            pen = tstep * 1e-6
            if i1p > i1pmax or i1p < i1pmin:
                pen += i1p**2
            if v1p > v1pmax or v1p < v1pmin:
                pen += v1pmax**2
            if v1next > v1max or v1next < v1min:
                pen += v1max**2
            return pen

        tstep = scipy.optimize.ftomin(tomin, 0)
        return tstep

    def set_vprop(self, vmin, vmax, vpmin, vpmax):
        self.v1max = vmax
        self.v1pmax = vpmax
        self.v1min = vmin
        self.v1pmin = vpmin


class LHCCircuits(object):
    def get_pc_names_from_lsa():
        lsa = pjlsa.LSAClient()
        pcs = lsa.findParameterNames(groupName="ALL MAGNETS", regexp=".*/IREF")
        pcs = [pc.split("/IREF")[0] for pc in pcs]
        return pcs

    def get_pc(self, pcs):
        out = {}
        pcs = set(pcs)
        for name in pcs:
            if name.endswith("B1") and name.replace("B1", "B2") in pcs:
                pc = Circuit2in1(name)
                out[pc.madname1] = pc
                out[pc.madname2] = pc
            elif name.endswith("B2"):
                pass
            else:
                pc = CircuitSingle(name)
                out[pc.madname] = pc
        return out

    @classmethod
    def from_dir(cls, dirname):
        pcs = {}
        for name in os.listdir(dirname):
            fn = os.path.join(dirname, name)
            if os.path.isfile(fn):
                pcs[name] = eval(open(fn).read(), {"nan": np.nan}, globals())
        return cls(pcs)

    def __init__(self, pcs):
        # self.pcs=self.get_pc(self.get_pc_names_from_lsa())
        self.pcs = pcs

    def store_directory(self, dirname):
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        for name, pc in self.pcs.items():
            fn = os.path.join(dirname, name)
            open(fn, "w").write(repr(pc))

    def load_directory(self, dirname):
        self.pcs = {}
        for name in os.listdir(dirname):
            fn = os.path.join(dirname, name)
            if os.path.isfile(fn):
                self.pcs[name] = eval(open(fn).read(), {"nan": np.nan}, globals())

    def set_calib_from_lsa(self):
        lsa = pjlsa.LSAClient()
        for name, pc in self.pcs.items():
            pc.set_calib_from_lsa(lsa=lsa)

    def set_pcinfo_from_lsa(self, fix=True):
        lsa = pjlsa.LSAClient()
        for name, pc in self.pcs.items():
            pc.set_pcinfo_from_lsa(lsa=lsa)
        if fix:
            self.set_vprop_from_pattern("RPMBA.*", -5, 5, -0.05, 0.05)  # trims
            self.set_vprop_from_pattern("RPMBB.*", -5, 2, -0.05, 0.05)  # ms
            self.set_vprop_from_pattern("RPHGA.*", 0.0, 5, -0.05, 0.05)  # q7-q10
            self.set_vprop_from_pattern("RPHGB.*", 0.0, 5, -0.05, 0.05)  # q5-q6
            self.set_vprop_from_pattern("RPHH.*", 0.0, 5, -0.05, 0.05)  # q4
            from pyoptics import madlang

            self.set_polarity_from_seq(madlang.open("lhc_as-built.seq"))
            self.scale_polarity_in_pattern("RPMBB.*RS[FD][12].*", 0.5)  # ms

    def set_pcinfo_from_pattern(self, regexp, attrname, value):
        creg = re.compile(regexp)
        for name, pc in self.pcs.items():
            if creg.match(name):
                print("set %s.%s=%s" % (name, attrname, value))
                setattr(pc, attrname, value)

    def set_vprop_from_pattern(self, regexp, vmin, vmax, vpmin, vpmax):
        creg = re.compile(regexp)
        for name, pc in self.pcs.items():
            if creg.match(pc.name):
                print("set vprop %-12s %s" % (name, pc.name))
                pc.set_vprop(vmin, vmax, vpmin, vpmax)

    def scale_polarity_in_pattern(self, regexp, scale):
        creg = re.compile(regexp)
        for name, pc in self.pcs.items():
            if creg.match(pc.name):
                print("scale polarity %-12s %s" % (name, pc.name))
                if hasattr(pc, "polarity1"):
                    pc.polarity1 *= scale
                if hasattr(pc, "polarity2"):
                    pc.polarity2 *= scale

    def mk_squeeze2(
        self,
        strtable,
        momentum,
        tstart,
        tstop,
        pcnames=None,
        tstepmin=0.02,
        tstepmax=10,
    ):
        # setup
        if pcnames is None:
            pcnames = set(self.pcs).intersection(set(strtable))
            pcnames = [
                p for p in pcnames if "b1" in p and p.replace("b1", "b2") in pcnames
            ]
        nsteps = len(strtable[pcnames[0]])
        # init
        tcurrent = 0
        tt = 0
        maxstep = 0.5
        tsteps = [tcurrent]
        isteps = {}
        for name in pcnames:
            pc = self.pcs[name]
            k1 = strtable[pc.madname1][0]
            k2 = strtable[pc.madname2][0]
            i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
            isteps[pc.madname1] = [i1]
            isteps[pc.madname2] = [i2]
        # loop at constant k
        while tcurrent < tstart - maxstep:
            tcurrent += maxstep
            tsteps.append(tcurrent)
            for name in pcnames:
                pc = self.pcs[name]
                k1 = strtable[pc.madname1][0]
                k2 = strtable[pc.madname2][0]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                # print(k1,k2,tsteps[-1],momentum(tsteps[-1]),i1,i2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
            # fmt="%4d: %7.3fs, %7.2f GeV"
            # print(fmt%(tt,tcurrent,momentum(tcurrent)/1e9))
        # loop varying k
        for name in pcnames:
            pc = self.pcs[name]
            isteps[pc.madname1].pop(-1)
            isteps[pc.madname2].pop(-1)
        for ttk in range(nsteps):
            # print(ttk,tcurrent)
            # find min step for all strengths
            for name in pcnames:
                maxstep = 1e-9
                pc = self.pcs[name]
                k1 = list(strtable[pc.madname1][ttk : ttk + 2])
                k2 = list(strtable[pc.madname2][ttk : ttk + 2])
                while len(k1) < 3:
                    k1 = k1 + [k1[-1]]
                    k2 = k2 + [k2[-1]]
                step = pc.get_min_tstep2(
                    momentum, tsteps[-2], k1, k2, maxstep, tstepmin, tstepmax
                )
                if step > maxstep:
                    maxstep = step
                    maxname = pc.madname1
                # tsteps_all[name].append(step)
            # choose step
            # steps=[(tsteps_all[name][-1],name) for name in  pcnames]
            # steps.sort(reverse=True)
            # print(steps)
            # maxstep,maxname=max(steps)
            # if ttk>nsteps-20:
            #    scale=1+(0.5*(ttk-nsteps+21)/9)**4
            #    print(ttk-nsteps+20, scale)
            #    maxstep*=scale
            # maxstep*=np.exp(float(ttk)/(nsteps-1))
            tsteps[-1] = tsteps[-2] + maxstep
            tsteps.append(tsteps[-1] + maxstep)
            tcurrent = tsteps[-2]
            # limits=' '.join(["%s:%.2fs"%(st[1],st[0]) for st in steps[:2]])
            fmt = "%4d: %7.3fs, %7.2f GeV, lim: %s:%5.3fs"
            print(fmt % (ttk, tcurrent, momentum(tcurrent) / 1e9, maxname, maxstep))
            # update step
            for name in pcnames:
                pc = self.pcs[name]
                k1 = strtable[pc.madname1][ttk]
                k2 = strtable[pc.madname2][ttk]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
        tsteps.pop(-1)
        for name in pcnames:
            pc = self.pcs[name]
            isteps[pc.madname1].pop(-1)
            isteps[pc.madname2].pop(-1)
        tsteps[-1] += 10
        # for name in pcnames:
        #   pc=self.pcs[name]
        #   k1=strtable[pc.madname1][ttk]
        #   k2=strtable[pc.madname2][ttk]
        #   i1,i2=pc.k2i(k1,k2,momentum(tcurrent))
        #   isteps[pc.madname1].append(i1)
        #   isteps[pc.madname2].append(i2)
        # loop at constant k
        # tsteps.pop(-1)
        # tsteps[-3]+=10
        # tsteps[-2]+=5
        # tsteps[-1]+=2
        for name in pcnames:
            pc = self.pcs[name]
            k1 = strtable[pc.madname1][ttk]
            k2 = strtable[pc.madname2][ttk]
            i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
            isteps[pc.madname1].append(i1)
            isteps[pc.madname2].append(i2)
        maxstep = 5.5
        print(tsteps[-1], k1)
        pp = True
        tcurrent = tsteps[-1]
        while tcurrent < tstop:
            tcurrent += maxstep
            if tcurrent > tstop:
                tcurrent = tstop
            tsteps.append(tcurrent)
            for name in pcnames:
                pc = self.pcs[name]
                k1 = strtable[pc.madname1][-1]
                if pp:
                    print("ccc", tcurrent, k1)
                    pp = False
                k2 = strtable[pc.madname2][-1]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                # print(k1,k2,tsteps[-1],momentum(tsteps[-1]),i1,i2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
            # fmt="%4d: %7.3fs, %7.2f GeV"
            # print(fmt%(tt,tcurrent,momentum(tcurrent)/1e9))
        print("name             Imin    Imax   |I'|max   |I''|max  Vmin    V'max")
        dsum = 0
        for name in sorted(pcnames):
            pc = self.pcs[name]
            pc.mk_simul(tsteps, isteps)
            dsum += pc.mk_simul_summ()
        print(f"Residual {dsum}")
        return tsteps, isteps

    def mk_squeeze(
        self,
        strtable,
        momentum,
        tstart,
        tstop,
        pcnames=None,
        tstepmin=0.02,
        w1=1,
        w2=0,
        w3=0,
    ):
        if pcnames is None:
            pcnames = set(self.pcs).intersection(set(strtable))
            pcnames = [
                p for p in pcnames if "b1" in p and p.replace("b1", "b2") in pcnames
            ]
        nsteps = len(strtable[pcnames[0]])
        # init
        tcurrent = 0
        tt = 0
        tsteps = [tcurrent]
        tsteps_all = {}
        isteps = {}
        ipsteps = {}
        for name in pcnames:
            pc = self.pcs[name]
            k1 = strtable[pc.madname1][0]
            k2 = strtable[pc.madname2][0]
            i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
            tsteps_all[name] = [0]
            isteps[pc.madname1] = [i1]
            isteps[pc.madname2] = [i2]
            ipsteps[pc.madname1] = [0]
            ipsteps[pc.madname2] = [0]
            # print(pc.madname1,k1,i1)
            # print(pc.madname2,k2,i2)
        # loop at constant k
        maxstep = 0.5
        while tcurrent < tstart - maxstep:
            tt += 1
            tcurrent += maxstep
            tsteps.append(tcurrent)
            for name in pcnames:
                pc = self.pcs[name]
                k1 = strtable[pc.madname1][0]
                k2 = strtable[pc.madname2][0]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                # print(k1,k2,tsteps[-1],momentum(tsteps[-1]),i1,i2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                i1p = (isteps[pc.madname1][tt] - isteps[pc.madname1][tt - 1]) / maxstep
                i2p = (isteps[pc.madname2][tt] - isteps[pc.madname2][tt - 1]) / maxstep
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
            # fmt="%4d: %7.3fs, %7.2f GeV"
            # print(fmt%(tt,tcurrent,momentum(tcurrent)/1e9))
        # loop varying k
        i1pforce = None
        i2pforce = None
        tt0 = tt
        for ttk in range(1, nsteps + 1):
            tt += 1
            # find min step for all strengths
            # print(tt-tt0-1)
            maxstep = 1e-9
            for name in pcnames:
                pc = self.pcs[name]
                if ttk < nsteps:
                    k1 = strtable[pc.madname1][ttk - 1]
                    k2 = strtable[pc.madname2][ttk - 1]
                    k1n = strtable[pc.madname1][ttk]
                    k2n = strtable[pc.madname2][ttk]
                else:
                    k1 = strtable[pc.madname1][-1]
                    k2 = strtable[pc.madname2][-1]
                    k1n = strtable[pc.madname1][-1]
                    k2n = strtable[pc.madname2][-1]
                i1p = ipsteps[pc.madname1][-1]
                i2p = ipsteps[pc.madname2][-1]
                step, i1p, i2p = pc.get_min_tstep(
                    momentum,
                    tcurrent,
                    k1,
                    k1n,
                    i1p,
                    k2,
                    k2n,
                    i2p,
                    i1pforce,
                    i2pforce,
                    maxstep,
                    tstepmin,
                )
                if step > maxstep:
                    maxstep = step
                    maxname = pc.madname1
                # tsteps_all[name].append(step)
            # choose step
            # steps=[(tsteps_all[name][-1],name) for name in  pcnames]
            # steps.sort(reverse=True)
            # print(steps)
            # maxstep,maxname=max(steps)
            # if ttk>nsteps-20:
            #    scale=1+(0.5*(ttk-nsteps+21)/9)**4
            #    print(ttk-nsteps+20, scale)
            #    maxstep*=scale
            # maxstep*=np.exp(float(ttk)/(nsteps-1))
            step1 = tsteps[-1] - tsteps[-2]
            step2 = tsteps[-2] - tsteps[-3]
            maxstep = step2 * w3 + step1 * w2 + maxstep * w1
            tcurrent += maxstep
            tsteps.append(tcurrent)
            # limits=' '.join(["%s:%.2fs"%(st[1],st[0]) for st in steps[:2]])
            fmt = "%4d: %7.3fs, %7.2f GeV, lim: %s:%5.3fs"
            print(fmt % (ttk, tcurrent, momentum(tcurrent) / 1e9, maxname, maxstep))
            # update step
            for name in pcnames:
                pc = self.pcs[name]
                # print(name,len(isteps[pc.madname1]),len(isteps[pc.madname2]))
                if ttk < nsteps:
                    k1 = strtable[pc.madname1][ttk]
                    k2 = strtable[pc.madname2][ttk]
                else:
                    k1 = strtable[pc.madname1][-1]
                    k2 = strtable[pc.madname2][-1]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                # print(k1,k2,tsteps[-1],momentum(tsteps[-1]),i1,i2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                i1p = (isteps[pc.madname1][-1] - isteps[pc.madname1][-2]) / maxstep
                i2p = (isteps[pc.madname2][-1] - isteps[pc.madname2][-2]) / maxstep
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
        # loop consant k
        # loop at constant k
        maxstep = 0.5
        while tcurrent < tstop:
            tt += 1
            tcurrent += maxstep
            if tcurrent > tstop:
                tcurrent = tstop
            tsteps.append(tcurrent)
            for name in pcnames:
                pc = self.pcs[name]
                k1 = strtable[pc.madname1][-1]
                k2 = strtable[pc.madname2][-1]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurrent))
                # print(k1,k2,tsteps[-1],momentum(tsteps[-1]),i1,i2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                i1p = (isteps[pc.madname1][-1] - isteps[pc.madname1][-2]) / maxstep
                i2p = (isteps[pc.madname2][-1] - isteps[pc.madname2][-2]) / maxstep
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
            # fmt="%4d: %7.3fs, %7.2f GeV"
            # print(fmt%(tt,tcurrent,momentum(tcurrent)/1e9))
        print("name             Imin    Imax   |I'|max   |I''|max  Vmin    V'max")
        dsum = 0
        for name in sorted(pcnames):
            pc = self.pcs[name]
            pc.mk_simul(tsteps, isteps)
            dsum += pc.mk_simul_summ()
        print(f"Residual {dsum}")
        return tsteps, isteps, ipsteps

    def mk_squeeze_tstep(
        self, strtable, momentum, tsteps, ksteps, pcnames=None, newt=None
    ):
        if pcnames is None:
            pcnames = set(self.pcs).intersection(set(strtable))
            pcnames = [
                p for p in pcnames if "b1" in p and p.replace("b1", "b2") in pcnames
            ]
        from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

        isteps = {}
        for name in pcnames:
            pc = self.pcs[name]
            isteps[pc.madname1] = []
            isteps[pc.madname2] = []
        for name in pcnames:
            pc = self.pcs[name]
            for tcurr, kcurr in zip(tsteps, ksteps):
                k1 = strtable[pc.madname1][kcurr]
                k2 = strtable[pc.madname2][kcurr]
                i1, i2 = pc.k2i(k1, k2, momentum(tcurr))
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
            if newt is not None:
                i1 = np.array(isteps[pc.madname1])
                i2 = np.array(isteps[pc.madname2])
                i1 = InterpolatedUnivariateSpline(tsteps, i1)(newt)
                i2 = InterpolatedUnivariateSpline(tsteps, i2)(newt)
                isteps[pc.madname1] = i1
                isteps[pc.madname2] = i2
        print("name             Imin    Imax   |I'|max   |I''|max  Vmin    V'max  summ")
        dsum = 0
        if newt is not None:
            tsteps = newt
        self.issues = {}
        for name in sorted(pcnames):
            pc = self.pcs[name]
            pc.mk_simul(tsteps, isteps)
            dsum += pc.mk_simul_summ()
            issues = pc.t[np.where(pc.mk_simul_summ_vec() > 0)[0]]
            if len(issues) > 0:
                self.issues[name] = issues
        print(f"Residual {dsum}")
        return tsteps, isteps, dsum

    def mk_squeeze_k(
        self,
        strtable,
        momentum,
        tstart,
        tstop,
        tstep,
        kstart=0.0,
        kstop=1.0,
        pcnames=None,
    ):
        if pcnames is None:
            pcnames = set(self.pcs).intersection(set(strtable))
            pcnames = [
                p for p in pcnames if "b1" in p and p.replace("b1", "b2") in pcnames
            ]
        # init
        nsteps = int(np.ceil((tstop - tstart) / tstep))
        tcurrent = 0
        kcurrent = kstart
        isteps = {}
        ipsteps = {}
        kfunc = {}
        ksteps_all = {}
        ksteps = [kcurrent]
        tsteps = [tcurrent]
        tt = 0
        for name in pcnames:
            pc = self.pcs[name]
            k1 = strtable.mk_function(pc.madname1)
            k2 = strtable.mk_function(pc.madname2)
            kfunc[pc.madname1] = k1
            kfunc[pc.madname2] = k2
            k1curr = k1(kcurrent)
            k2curr = k2(kcurrent)
            pccurr = momentum(tcurrent)
            i1, i2 = pc.k2i(k1curr, k2curr, momentum(tcurrent))
            isteps[pc.madname1] = [i1]
            isteps[pc.madname2] = [i2]
            ipsteps[pc.madname1] = [0]
            ipsteps[pc.madname2] = [0]
            ksteps_all[name] = [0]

        # def update step
        def update_step(tnext, knext, pcnext):
            for name in pcnames:
                pc = self.pcs[name]
                k1 = kfunc[pc.madname1]
                k2 = kfunc[pc.madname2]
                i1, i2 = pc.k2i(k1(knext), k2(knext), pcnext)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                i1p = (isteps[pc.madname1][-1] - isteps[pc.madname1][-2]) / tstep
                i2p = (isteps[pc.madname2][-1] - isteps[pc.madname2][-2]) / tstep
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
                # print(pc.madname1,i1p,i2p)
            tsteps.append(tnext)
            ksteps.append(knext)

        # loop at constant k
        while tcurrent < tstart - tstep:
            tnext = tcurrent + tstep
            pcnext = momentum(tnext)
            update_step(tnext, kcurrent, pcnext)
            tcurrent = tnext
        fmt = "%4d: %7.3fs, %7.2f GeV, %7.3f"
        print(fmt % (len(tsteps), tnext, pcnext / 1e9, kcurrent))
        # print(len(tsteps),len(isteps[pc.madname1]))
        i1pforce = None
        i2pforce = None

        # loop varying k
        while kcurrent <= kstop and tcurrent < tstop:
            # find min k step for all strengths
            minknext = kstop
            tnext = tcurrent + tstep
            pccurrent = momentum(tcurrent)
            pcnext = momentum(tnext)
            for name in pcnames:
                pc = self.pcs[name]
                k1 = kfunc[pc.madname1]
                k2 = kfunc[pc.madname2]
                i1p = ipsteps[pc.madname1][-1]
                i2p = ipsteps[pc.madname2][-1]
                knext, i1p, i2p = pc.get_min_kstep(
                    pccurrent,
                    pcnext,
                    k1,
                    k2,
                    kcurrent,
                    minknext,
                    tstep,
                    i1p,
                    i2p,
                    i1pforce,
                    i2pforce,
                )
                if knext < minknext:
                    minknext = knext
                    minname = name
            knext = minknext
            fmt = "%4d: %7.3fs, %7.2f GeV, %7.3f: lim %s"
            print(fmt % (tt, tnext, pcnext / 1e9, knext, minname))
            update_step(tnext, knext, pcnext)
            tcurrent = tnext
            kcurrent = knext

        # loop at constant k
        while tcurrent < tstop:
            tnext = tcurrent + tstep
            if tnext > tstop:
                tnext = tstop
            pcnext = momentum(tnext)
            update_step(tnext, knext, pcnext)
            tcurrent = tnext

        # check results
        for name in pcnames:
            pc = self.pcs[name]
            print(name)
            pc.mk_simul(tsteps, isteps)
        return tsteps, ksteps, isteps, ipsteps

    def mk_squeeze_cubic(
        self, strtable, momentum, tstart, tstop, pcnames=None, tstepmax=1.0
    ):
        if pcnames is None:
            pcnames = set(self.pcs).intersection(set(strtable))
            pcnames = [
                p for p in pcnames if "b1" in p and p.replace("b1", "b2") in pcnames
            ]
        nsteps = len(strtable[pcnames[0]])
        # init
        cubics1 = {}
        cubics2 = {}
        isteps = {}
        ipsteps = {}
        tsteps = []
        tcurrent = 0
        i1p = 0
        i2p = 0
        ik = 0
        for name in pcnames:
            pc = self.pcs[name]
            cubics1[name] = []
            cubics2[name] = []
            isteps[pc.madname1] = []
            isteps[pc.madname2] = []
            ipsteps[pc.madname1] = []
            ipsteps[pc.madname2] = []
        tstep = tstepmax
        while tcurrent < tstart:
            for name in pcnames:
                pc = self.pcs[name]
                i1v = []
                i2v = []
                k1 = strtable[pc.madname1][0]
                k2 = strtable[pc.madname2][0]
                for ii in [0, 1, 2]:
                    i1, i2 = pc.k2i(k1, k2, momentum(tcurrent + tstep * ii))
                    i1v.append(i1)
                    i2v.append(i2)
                c1 = cubic_fit2(i1v[0], i1v[1], i1v[2], i1p, tstep, 2 * tstep)
                c2 = cubic_fit2(i2v[0], i2v[1], i2v[2], i2p, tstep, 2 * tstep)
                cubics1[name].append((tstep,) + c1)
                cubics2[name].append((tstep,) + c2)
                i1, i1p, i1pp = cubic_val(0, *c1)
                i2, i2p, i2pp = cubic_val(0, *c2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
                i1, i1p, i1pp = cubic_val(tstep, *c1)
                i2, i2p, i2pp = cubic_val(tstep, *c2)
            pccurrent = momentum(tcurrent)
            fmt = "%4d: %7.3fs, %7.2f GeV"
            print(fmt % (len(tsteps), tcurrent, pccurrent / 1e9))
            tsteps.append(tcurrent)
            tcurrent += tstep
        # loop at varying k
        for ttk in range(0, nsteps):
            # optimize tstep
            tstepmin = 0.2
            for name in pcnames:
                pc = self.pcs[name]
                pen = 1
                tstep = 0.0
                while pen > 0 and tstep < 30:
                    tstep += 0.2
                    i1v = []
                    i2v = []
                    for ii in [0, 1, 2]:
                        ttkii = ttk + ii
                        if ttkii >= nsteps:
                            ttkii = nsteps - 1
                        k1 = strtable[pc.madname1][ttkii]
                        k2 = strtable[pc.madname2][ttkii]
                        i1, i2 = pc.k2i(k1, k2, momentum(tcurrent + tstep * ii))
                        i1v.append(i1)
                        i2v.append(i2)
                    c1 = cubic_fit2(i1v[0], i1v[1], i1v[2], i1p, tstep, 2 * tstep)
                    c2 = cubic_fit2(i2v[0], i2v[1], i2v[2], i2p, tstep, 2 * tstep)
                    pen = pc.model_get_constraints_cubic(c1, c2, tstep)
                # print("%s %s %s"%(name,tstep,tstepmin))
                if tstep > tstepmin:
                    tstepmin = tstep
                    namemin = name
            # calculate  mintstep
            # chosen tstep
            for name in pcnames:
                pc = self.pcs[name]
                i1v = []
                i2v = []
                for ii in [0, 1, 2]:
                    ttkii = ttk + ii
                    if ttkii >= nsteps:
                        ttkii = nsteps - 1
                    k1 = strtable[pc.madname1][ttkii]
                    k2 = strtable[pc.madname2][ttkii]
                    i1, i2 = pc.k2i(k1, k2, momentum(tcurrent + tstep * ii))
                    i1v.append(i1)
                    i2v.append(i2)
                # TODO:optimize tstep
                c1 = cubic_fit2(i1v[0], i1v[1], i1v[2], i1p, tstep, 2 * tstep)
                c2 = cubic_fit2(i2v[0], i2v[1], i2v[2], i2p, tstep, 2 * tstep)
                cubics1[name].append((tcurrent, tstep) + c1)
                cubics2[name].append((tcurrent, tstep) + c2)
                i1, i1p, i1pp = cubic_val(0, *c1)
                i2, i2p, i2pp = cubic_val(0, *c2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
                i1, i1p, i1pp = cubic_val(tstep, *c1)
                i2, i2p, i2pp = cubic_val(tstep, *c2)
            pccurrent = momentum(tcurrent)
            fmt = "%4d: %7.3fs, %7.2f GeV, %7.3f: lim %s"
            minname = "none"
            print(fmt % (len(tsteps), tcurrent, pccurrent / 1e9, tstep, namemin))
            tsteps.append(tcurrent)
            tcurrent += tstep
        # loop at constant k
        tstep = tstepmax
        while tcurrent <= tstop + tstep:
            for name in pcnames:
                pc = self.pcs[name]
                i1v = []
                i2v = []
                k1 = strtable[pc.madname1][-1]
                k2 = strtable[pc.madname2][-1]
                for ii in [0, 1, 2]:
                    i1, i2 = pc.k2i(k1, k2, momentum(tcurrent + tstep * ii))
                    i1v.append(i1)
                    i2v.append(i2)
                c1 = cubic_fit2(i1v[0], i1v[1], i1v[2], i1p, tstep, 2 * tstep)
                c2 = cubic_fit2(i2v[0], i2v[1], i2v[2], i2p, tstep, 2 * tstep)
                cubics1[name].append((tcurrent, tstep) + c1)
                cubics2[name].append((tcurrent, tstep) + c2)
                i1, i1p, i1pp = cubic_val(0, *c1)
                i2, i2p, i2pp = cubic_val(0, *c2)
                isteps[pc.madname1].append(i1)
                isteps[pc.madname2].append(i2)
                ipsteps[pc.madname1].append(i1p)
                ipsteps[pc.madname2].append(i2p)
                i1, i1p, i1pp = cubic_val(tstep, *c1)
                i2, i2p, i2pp = cubic_val(tstep, *c2)
            pccurrent = momentum(tcurrent)
            fmt = "%4d: %7.3fs, %7.2f GeV"
            print(fmt % (len(tsteps), tcurrent, pccurrent / 1e9))
            tsteps.append(tcurrent)
            if tcurrent == tstop:
                break
            tcurrent += tstep
            if tcurrent > tstop:
                tstep = tcurrent - tstop
                tcurrent = tstop

        # post processing
        print("name             Imin    Imax   |I'|max   Vmin    V'max   I''max")
        dsum = 0
        for name in sorted(pcnames):
            pc = self.pcs[name]
            pc.mk_simul(tsteps, isteps)
            dsum += pc.mk_simul_summ()
        print(f"Residual {dsum}")
        return tsteps, isteps, ipsteps

    def mk_simul(self, tsteps, isteps):
        for name, pc in self.pcs.items():
            if name in isteps:
                pc.mk_simul(tsteps, isteps)

    def set_polarity_from_seq(self, seq):
        deps = seq.build_dep()
        for name, pc in self.pcs.items():
            if hasattr(pc, "name2"):
                name1 = name.replace(".", "_")
                name2 = name1.replace("b1", "b2")
                if name1 in deps:
                    pc.polarity1 = seq[deps[name1][0][0]].polarity
                if name2 in deps:
                    pc.polarity2 = seq[deps[name2][0][0]].polarity
            else:
                name1 = name.replace(".", "_")
                if name1 in deps:
                    pc.polarity1 = seq[deps[name1][0][0]].polarity

    def mk_hllhc(self):
        pcs = self.pcs
        for irn in "15":
            for lr in "lr":
                q5 = "kq5.%s%sb1" % (lr, irn)
                q4 = "kq4.%s%sb1" % (lr, irn)
                pc = pcs[q4].copy()
                pc.madname1 = q5
                pc.madname2 = q5.replace("b1", "b2")
                pc.polarity1 *= -1
                pc.polarity2 *= -1
                pcs[q5] = pc
        return self.__class__(pcs)
