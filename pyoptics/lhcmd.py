import re

import pytimber
from . import rdmdate
from .sddsdata import sddsdata

from numpy import *
import numpy as np
from matplotlib.pyplot import *

from .picker import Picker


def gen_both(s):
    regexp = re.compile(s, re.IGNORECASE)
    out = []
    for l in open("powerlist"):
        res = regexp.match(l.strip())
        if res:
            name = res.group()
            out.append(name)
            out.append(name.replace("I_MEAS", "I_REF"))
    return out


def gen_iref(s):
    regexp = re.compile(s, re.IGNORECASE)
    out = []
    for l in open("powerlist"):
        res = regexp.match(l.strip())
        if res:
            name = res.group()
            out.append(name.replace("I_MEAS", "I_REF"))
    return out


def gen_imeas(s):
    regexp = re.compile(s, re.IGNORECASE)
    out = []
    for l in open("powerlist"):
        res = regexp.match(l.strip())
        if res:
            name = res.group()
            out.append(name)
    return out


def getdata(vs, t1="2011-05-04 14:38:00.000", t2="2011-05-04 14:45:00.000"):
    ldb = pytimber.LoggingDB()
    data = ldb.get(",".join(vs), t1, t2, conf="mdb.conf")
    # print '\n'.join(data['log'])
    for k in data["datavars"]:
        data[k][0] /= 1000.0
    return data


def repknobs(data, dateint, vs):
    for name in vs:
        if "I_REF" in name:
            t, v = data[name]
            tn = rdmdate.parsedate(dateint[0])
            print("%-24s" % name.split(":")[0], end=" ")
            try:
                tref = where(t > tn)[0][0]
                print("%10.6f" % v[tref], end=" ")
                for ti in dateint[1:]:
                    tn = rdmdate.parsedate(ti)
                    ti = where(t > tn)[0][0]
                    print("%11.6f" % (v[ti] - v[tref]), end=" ")
                print()
            except:
                print("missing")


def allrepknobs(t1, t2, dateint, *vsl):
    print("\n".join(dateint))
    for vs in vsl:
        data = getdata(vs, t1=t1, t2=t2)
        repknobs(data, dateint, vs)


def depick(p):
    print(
        "dateint=%s"
        % repr([rdmdate.dumpdate(ti, fmt="%Y-%m-%d %H:%M:%S") for ti, vi in p.data])
    )


# arcs
mb = gen_both(".*RB\.A.*")
mcs = gen_both(".*RCS\.A.*")
mco = gen_both(".*RCO\.A.*")
mcd = gen_both(".*RCD\.A.*")
mq = gen_both(".*RQ[FD].*")

mall = gen_both(".*")

# knobs tune chroma coupling
mqt = gen_both(".*RQT[FD].*")
mqtfb1 = gen_both(".*RQTF.*B1.*")
mqtdb1 = gen_both(".*RQTD.*B1.*")
mqtfb2 = gen_both(".*RQTF.*B2.*")
mqtdb2 = gen_both(".*RQTD.*B2.*")
ms = gen_both(".*RS[FD].*")
msfb1 = gen_both(".*RSF.*B1.*")
msdb1 = gen_both(".*RSD.*B1.*")
msfb2 = gen_both(".*RSF.*B2.*")
msdb2 = gen_both(".*RSD.*B2.*")
mqs = gen_both(".*RQS\..*")
mqsb1 = gen_both(".*RQS\..*B1.*")
mqsb2 = gen_both(".*RQS\..*B2.*")

# crossing scheme
mcbx_ir1 = gen_both(".*RCBX.*[LR]1.*")
mcbx_ir2 = gen_both(".*RCBX.*[LR]2.*")
mcbx_ir5 = gen_both(".*RCBX.*[LR]5.*")
mcbx_ir8 = gen_both(".*RCBX.*[LR]8.*")


# squeeze
ir1 = gen_both(".*RQ.*R1B[12].*") + gen_both(".*RT?QX.*[LR]1.*")
ir2 = gen_both(".*RQ.*R2B[12].*") + gen_both(".*RT?QX.*[LR]2.*")
ir5 = gen_both(".*RQ.*R5B[12].*") + gen_both(".*RT?QX.*[LR]5.*")
ir8 = gen_both(".*RQ.*R8B[12].*") + gen_both(".*RT?QX.*[LR]8.*")


def orbk(s):
    res = re.search(r"([4-6])\.([LR]).*", s)
    n, l = res.groups()
    if l == "R":
        rn = int(n) + 500
    else:
        rn = int(n) + 200
    # print s,rn
    return rn


mcb_ir1 = sorted(gen_both(".*RR.*RCB.*[4-6]\.[LR]1.*"), key=orbk)
mcb_ir2 = sorted(gen_both(".*UA.*RCB.*[4-6]\.[LR]2.*"), key=orbk)
mcb_ir5 = sorted(gen_both(".*RR.*RCB.*[4-6]\.[LR]5.*"), key=orbk)
mcb_ir8 = sorted(gen_both(".*UA.*RCB.*[4-6]\.[LR]8.*"), key=orbk)


bctdc = [
    "LHC.BCTDC.A6R4.B1:BEAM_INTENSITY",
    "LHC.BCTDC.B6R4.B1:BEAM_INTENSITY",
    "LHC.BCTDC.A6R4.B2:BEAM_INTENSITY",
    "LHC.BCTDC.B6R4.B2:BEAM_INTENSITY",
]

bctfr = [
    "LHC.BCTFR.A6R4.B1:BEAM_INTENSITY",
    "LHC.BCTFR.A6R4.B2:BEAM_INTENSITY",
]

rffreq = ["ALB.SR4.B1:FGC_FREQ", "ALB.SR4.B2:FGC_FREQ"]


def mkdate(t):
    return rdmdate.dumpdate(t, fmt="%Y-%m-%d\n%H:%M:%S.SSS")


def myplot(data):
    clf()
    for vi in range(0, len(data["datavars"]), 2):
        name = data["datavars"][vi]
        nameref = data["datavars"][vi + 1]
        t, v = data[vi]
        plot(t, v, ".-", label=name)
        t, v = data[vi + 1]
        plot(t, v, "-", label=nameref)
    myticks()


def myticks():
    t1, t2 = xlim()
    dt = linspace(t1, t2, 10)
    lt = [mkdate(tt) for tt in dt]
    xticks(dt, lt)


import matplotlib.dates
import calendar


def bpmdata_date(fn):
    p = fn.split(".")[0].split("@")
    year, month, day = list(map(int, p[2].split("_")))
    hour, minute, second, msecond = list(map(int, p[3].split("_")))
    timestamp = calendar.timegm((year, month, day, hour, minute, second, 0, -1))
    timestamp += msecond * 1e-3
    return timestamp


def repdiff(ref, data):
    for vs in ref["datavars"]:
        if vs.endswith("I_REF"):
            rt, rv = ref[vs]
            dt, dv = data[vs]
            rv = rv[0]
            dv = dv[0]
            if rv != dv:
                t1 = rdmdate.dumpdate(rt[0])
                t2 = rdmdate.dumpdate(dt[0])
                print("%-30s %10g %10g %10s %10s" % (vs, rv, rv - dv, t1, t2))
    print()


(
    b1h,
    b1v,
    b2h,
    b2v,
) = "LHC.BQBBQ.UA47.FFT1_B1:FFT_DATA_H LHC.BQBBQ.UA47.FFT1_B1:FFT_DATA_V LHC.BQBBQ.UA43.FFT1_B2:FFT_DATA_H LHC.BQBBQ.UA43.FFT1_B2:FFT_DATA_V".split()


frfb1 = "ALB.SR4.B1:FGC_FREQ"
frfb2 = "ALB.SR4.B2:FGC_FREQ"

import gzip
import os
from tempfile import mktemp

from scipy import optimize
from numpy.fft import *
import matplotlib.pyplot as pl

from .pydspro import t2f
from . import harmonic_fit


class betabeatBPM(object):
    def __init__(self, fn):
        self.filename = fn
        fh = gzip.open(fn)
        self.bunchid = int(fh.readline().split(":")[1])
        self.turns = int(fh.readline().split(":")[1])
        self.monitors = int(fh.readline().split(":")[1])
        flags = []
        bpms = []
        spos = []
        data = []
        for line in fh:
            line = line.split()
            flags.append(int(line[0]))
            bpms.append(line[1])
            spos.append(line[2])
            data.append(line[3:])
            # assert len(data[-1])==self.turns
        self.flags = array(flags, dtype=bool)
        self.spos = array(spos, dtype=float)
        self.data = array(data, dtype=float)
        self.bpms = array(bpms)
        self.bad = self.find_bad_bpms()
        self.t = arange(self.turns)

    def find_bad_bpms(self):
        bad = self.find_bpm_with_null_values()
        bad |= self.find_bpm_constant_values()
        bad |= self.find_bpm_huge_values()
        return bad

    def find_bpm_with_null_values(self, n=200):
        cnd = array([len(where(v == 0)[0]) > n for v in self.data])
        print("BPM:%d bpms have >%d values = 0" % (sum(cnd), n))
        return cnd

    def find_bpm_constant_values(self, n=200):
        cnd = array([len(where(diff(v) == 0)[0]) > 200 for v in self.data])
        print("BPM:%d bpms have >%d with same values" % (sum(cnd), n))
        return cnd

    def find_bpm_huge_values(self, n=1000):
        cnd = any(abs(self.data) > 10, axis=1)
        print("BPM:%d bpms have >%g values" % (sum(cnd), n))
        return cnd

    def plot_bad_bpms(m):
        out = []
        for i in where(m.bad)[0]:
            pl.clf()
            pl.plot(m.data[i])
            if m.xidx[i]:
                pl.title("%s x i=%d" % (m.bpms[i], i))
            else:
                pl.title("%s y i=%d" % (m.bpms[i], i))
            fn = "bad_%d.pdf" % i
            pl.savefig(fn)
            out.append(fn)
        fnames = " ".join(out)
        os.system("pdftk %s cat output %s_bad.pdf" % (fnames, m.title()))
        os.system("rm %s" % fnames)

    def title(self):
        return os.path.split(self.filename)[1].split(".sdds.gz")[0]

    def s2bpm(self, s):
        return self.bpms[self.spos == s]

    def mk_fitlsq(self):
        t = self.t
        out = []
        for i in range(len(self.data)):
            if not self.bad[i]:
                co, ff, a, p, res = lsqmax(self.data[i], t)
                s = self.spos[i]
                out.append([s, ff, a, p, co, res, i, self.bpms[i]])
        self.fit_tbl = betabeatFit(*list(zip(*sorted(out))))

    def refine_fit(self, qx=0.27, qy=0.32, qtol=0.02):
        self.fitx = self.fit_tbl.filter_tune(qx - qtol / 2, qx + qtol / 2)
        self.fity = self.fit_tbl.filter_tune(qy - qtol / 2, qy + qtol / 2)

    def plot_bety(m, t1):
        m._plot_bet(m.fity, t1, t1.bety, "y")

    def plot_betx(m, t1):
        m._plot_bet(m.fitx, t1, t1.betx, "x")

    def _plot_bet(m, fit, t1, betA, lbl):
        idx = array([where(t1.name == bp)[0][0] for bp in fit.bpms])
        betA = betA[idx]
        betB = fit.amp**2
        fact = (sqrt(betA / betB)).mean() ** 2
        # pl.figure()
        pl.plot(fit.s, betA, "b-", label=r"model $\beta_%s$" % lbl)
        pl.plot(fit.s, betB * fact, "r-", label=r"meas. $%g \cdot J_%s$" % (fact, lbl))
        pl.legend()
        pl.grid(True)
        pl.title(os.path.split(m.filename)[1])

    def plot_muy(m, t1):
        m._plot_mu(m.fity, t1, t1.muy, "y")

    def plot_mux(m, t1):
        m._plot_mu(m.fitx, t1, t1.mux, "x")

    def _plot_mu(m, fit, t1, muA, lbl):
        idx = array([where(t1.name == bp)[0][0] for bp in fit.bpms])
        muA = rad2deg(diff(unwrap(muA[idx] * 2 * pi)))
        muB = rad2deg(diff(unwrap(fit.phase)))
        pl.figure()
        pl.plot(fit.s[:-1], muA, "b-", label=r"model $\Delta \mu_%s$" % lbl)
        pl.plot(fit.s[:-1], muB, "r-", label=r"meas. $\Delta \mu_%s$" % lbl)
        pl.legend()
        pl.grid(True)
        pl.title(os.path.split(m.filename)[1])


class betabeatFit(object):
    def __init__(self, s, ff, a, p, co, res, iii, bpm):
        self.s = array(s)
        self.tune = array(ff)
        self.amp = array(a)
        self.phase = array(p)
        self.co = array(co)
        self.res = array(res)
        self.iii = array(iii)
        self.bpms = array(bpm)

    def filter_tune(self, qa, qb):
        idx = (self.tune <= qb) & (self.tune >= qa)
        tuneavg = self.tune[idx].mean()
        qmin = self.tune[idx].min()
        qmax = self.tune[idx].max()
        print(
            "Found %d bpm with tune %9.6f in (%9.6f,%9.6f)"
            % (sum(idx), tuneavg, qmin, qmax)
        )
        res = betabeatFit.__new__(betabeatFit)
        for k, v in list(self.__dict__.items()):
            setattr(res, k, v[idx])
        return res


def ham(t, ff, a, p):
    return a * cos(2 * pi * ff * t + p)


def getfftmax(v, f):
    fv = rfft(v)
    fv[0] = 0
    idx = abs(fv).argmax()
    ff = f[idx]
    a = abs(fv[idx]) / (len(v) / 2)
    p = angle(fv[idx])
    return ff, a, p


def lsqmax(v, t):
    f = t2f(t)
    x = getfftmax(v, f)
    co = v.mean()

    def ftomin(x):
        ff, a, p = x
        return sum((v - ham(t, ff, a, p) - co) ** 2)

    ff, a, p = optimize.fmin(ftomin, x, disp=True, xtol=1e-9, ftol=1e-9)
    res = ftomin([ff, a, p])
    return co, ff, a, p, sqrt(res)


def lsqmax2(v1, v2, t):
    f = t2f(t)
    n1 = float(len(v1))
    n2 = float(len(v2))
    # fftmax 1
    fv1 = rfft(v1)
    co1 = fv1[0] / n1
    idx1 = abs(fv1).argmax()
    f1 = idx1 / n1
    a11 = abs(fv1[idx1]) * 2 / n1
    p11 = angle(fv1[idx1])
    # fftmax 2
    fv2 = rfft(v2)
    co2 = fv2[0] / n2
    idx2 = abs(fv2).argmax()
    f2 = f[idx2] / n2
    a22 = abs(fv2[idx2]) * 2 / n2
    p22 = angle(fv2[idx2])
    # cross terms
    a12 = abs(fv1[idx2]) * 2 / n1
    p12 = angle(fv1[idx2])
    a21 = abs(fv2[idx1]) * 2 / n2
    p21 = angle(fv2[idx1])
    x = f1, f2, a11, a12, a21, a22, p11, p12, p21, p22

    def ftomin(x):
        f1, f2, a11, a12, a21, a22, p11, p12, p21, p22 = x
        t1 = 2 * pi * f1 * t
        t2 = 2 * pi * f2 * t
        vv1 = co1 + a11 * cos(t1 + p11) + a12 * cos(t2 + p12)
        vv2 = co2 + a21 * cos(t1 + p21) + a22 * cos(t2 + p22)
        return sum((vv1 - v1) ** 2 + (vv2 - v2) ** 2)

    x = optimize.fmin(ftomin, x, disp=True, xtol=1e-9, ftol=1e-9)
    f1, f2, a11, a12, a21, a22, p11, p12, p21, p22 = x
    res = ftomin(x)
    return co1, co2, f1, f2, a11, a12, a21, a22, p11, p12, p21, p22, res


class LHCBPM(object):
    def __repr__(self):
        if hasattr(self, "name"):
            return "<LHCBPM '%s'>" % (self.name)
        else:
            return object.__repr__(self)

    def __init__(self, fn):
        sdds = sddsdata(fn)
        self.filename = fn
        self.bpms = r_[sdds.data[0]["bpmNames"], sdds.data[0]["bpmNames"]]
        self.turns = sdds.data[0]["nbOfCapTurns"][0]
        xdata = sdds.data[0]["horPositionsConcentratedAndSorted"]
        ydata = sdds.data[0]["verPositionsConcentratedAndSorted"]
        self.data = r_[xdata, ydata].reshape(len(self.bpms), self.turns)
        self.xidx = zeros(len(self.bpms), dtype=bool)
        self.xidx[: len(self.bpms) / 2] = True
        self.yidx = ~self.xidx
        self.bad = self.find_bad_bpms()
        self.badxy = self.mk_badxy()
        self.t = arange(self.turns)
        self.goodx = (~self.badxy) & self.xidx
        self.goody = (~self.badxy) & self.yidx

    def find_bad_bpms(self):
        bad = self.find_bpm_with_null_values()
        bad |= self.find_bpm_constant_values()
        bad |= self.find_bpm_huge_values()
        return bad

    def mk_badxy(self):
        """
        depends: bad
        provides: badxy
        """
        if not hasattr(self, "badxy"):
            badxy = zeros(len(self.bpms), dtype=bool)
            for n in set(self.bpms[self.bad]):
                badxy[self.bpms == n] = True
            self.badxy = badxy
        return self.badxy

    def find_bpm_with_null_values(self, n=200):
        cnd = array([len(where(v == 0)[0]) > n for v in self.data])
        print("BPM:%d bpms have >%d values = 0" % (sum(cnd), n))
        return cnd

    def find_bpm_constant_values(self, n=200):
        cnd = array([len(where(diff(v) == 0)[0]) > 200 for v in self.data])
        print("BPM:%d bpms have >%d with same values" % (sum(cnd), n))
        return cnd

    def find_bpm_huge_values(self, n=30):
        cnd = any(abs(self.data) > 10, axis=1)
        print("BPM:%d bpms have >%g values" % (sum(cnd), n))
        return cnd

    def plot_bad_bpms(m):
        out = []
        for i in where(m.bad)[0]:
            pl.clf()
            pl.plot(m.data[i])
            pl.title("%s  i=%d" % (m.bpms[i], i))
            fn = "bad_%d.pdf" % i
            pl.savefig(fn)
            out.append(fn)
        fnames = " ".join(out)
        os.system("pdftk %s cat output %s_bad.pdf" % (fnames, m.title()))
        os.system("rm %s" % fnames)

    def title(self):
        return os.path.split(self.filename)[1].split(".sdds.gz")[0]

    def s2bpm(self, s):
        return self.bpms[self.spos == s]

    def mk_spos(self, twisstable):
        """
        provide: spos betx bety mux muy dx xidxs yidxs
        """
        self._mkfloatvec("spos betx bety mux muy dx")
        sign = 1
        if twisstable.header["sequence"] == "LHCB2":
            sign = -1
        for i, n in enumerate(self.bpms):
            try:
                self.spos[i] = twisstable.s[twisstable.name == n]
                self.betx[i] = twisstable.betx[twisstable.name == n]
                self.bety[i] = twisstable.bety[twisstable.name == n]
                self.dx[i] = twisstable.dx[twisstable.name == n]
                self.mux[i] = sign * twisstable.mux[twisstable.name == n]
                self.muy[i] = sign * twisstable.muy[twisstable.name == n]
            except ValueError:
                print("Error on %s in line %d" % (n, i))
        idx = where(self.xidx & ~self.badxy)[0]
        svec, sidx = list(zip(*sorted([(self.spos[i], i) for i in idx])))
        self.xidxs = array(sidx)
        idx = where(self.yidx & ~self.badxy)[0]
        svec, sidx = list(zip(*sorted([(self.spos[i], i) for i in idx])))
        self.yidxs = array(sidx)

    def _mkfloatvec(self, namelst):
        for n in namelst.split():
            setattr(self, n, np.zeros(len(self.data), dtype=float))

    def _setrow(self, namelst, i, *args):
        for n, v in zip(namelst.split(), args):
            getattr(self, n)[i] = v

    def mk_fitlsq(self):
        """
        depends: bad
        provides: tune amp phase co res
        """
        self._mkfloatvec("tune amp phase co res")
        for i in range(len(self.data)):
            if not self.bad[i]:
                print(self.bpms[i])
                co, ff, a, p, res = fit_single_lsq(self.data[i])
                self._setrow("co tune amp phase res", i, co, ff, a, p, res)

    def mk_coupled_coeff(self):
        """
        depends: bad
        provides: tune amp phase co res
        """
        t = self.t
        self._mkfloatvec("amp2 phase2")
        tunex = mean(self.tune[self.goodx])
        tuney = mean(self.tune[self.goody])
        for i in self.xidxs:
            a12, p12 = get_harm_coef(self.data[i], tuney)
            self._setrow("amp2 phase2", i, a12, p12)
        for i in self.yidxs:
            a21, p21 = get_harm_coef(self.data[i], tunex)
            self._setrow("amp2 phase2", i, a21, p21)

    def mk_fitlsq2(self):
        """
        depends: bad
        provides: tune amp phase co res
        """
        t = self.t
        vecname = "co tune amp phase co amp2 phase2 res"
        self._mkfloatvec(vecname)
        xybpms = len(self.data) / 2
        for i in range(xybpms):
            if not self.badxy[i]:
                print("fitting %s" % self.bpms[i])
                v1 = self.data[i]
                v2 = self.data[i + xybpms]
                x = fit_coupled_lsq(v1, v2)
                co1, co2, f1, f2, a11, a12, a21, a22, p11, p12, p21, p22, res = x
                self._setrow(vecname, i, co1, f1, a11, p11, a12, p12, res)
                self._setrow(vecname, i + xybpms, co2, f2, a22, p22, a21, p21, res)

    def plot_2dtunes(m, samescale=True, errorbar=None):
        tunex = m.tune[m.xidx & ~m.badxy]
        tuney = m.tune[m.yidx & ~m.badxy]
        pl.plot(tunex, tuney, ".")
        if errorbar is not None:
            errx = m.res[m.xidx & ~m.badxy] * errorbar
            erry = m.res[m.xidx & ~m.badxy] * errorbar
            pl.errorbar(tunex, tuney, xerr=errx, yerr=erry, fmt=".")
        if samescale:
            avgx = (tunex.max() + tunex.min()) / 2
            avgy = (tuney.max() + tuney.min()) / 2
            spdx = (tunex.max() - tunex.min()) / 2
            spdy = (tuney.max() - tuney.min()) / 2
            spd = max(spdx, spdy)
            pl.xlim(avgx - spd, avgx + spd)
            pl.ylim(avgy - spd, avgy + spd)
        return tunex, tuney

    def plot_betx(m):
        betA = sqrt(m.betx[m.xidxs])
        betB = m.amp[m.xidxs]
        k = mean(betA / betB) ** 2
        pl.plot(m.spos[m.xidxs], betA**2, "-", label=r"model $\beta_x$")
        pl.plot(m.spos[m.xidxs], betB**2 * k, "-", label=r"meas $\beta_x$")
        xlabel("s [m]")
        ylabel(r"$\beta$ [m]")

    def plot_bety(m):
        betA = sqrt(m.bety[m.yidxs])
        betB = m.amp[m.yidxs]
        k = mean(betA / betB) ** 2
        pl.plot(m.spos[m.yidxs], betA**2, "-", label=r"model $\beta_y$")
        pl.plot(m.spos[m.yidxs], betB**2 * k, "-", label=r"meas $\beta_y$")
        xlabel("s [m]")
        ylabel(r"$\beta$ [m]")

    def plot_betx_beat(m):
        betA = sqrt(m.betx[m.xidxs])
        betB = m.amp[m.xidxs]
        k = mean(betA / betB) ** 2
        pl.plot(
            m.spos[m.xidxs],
            betB**2 * k / betA**2 - 1,
            "-",
            label=r"model $\beta_x$",
        )
        xlabel("s [m]")
        ylabel(r"$\Delta\beta$")

    def plot_bety_beat(m):
        betA = sqrt(m.bety[m.yidxs])
        betB = m.amp[m.yidxs]
        k = mean(betA / betB) ** 2
        pl.plot(
            m.spos[m.yidxs],
            betB**2 * k / betA**2 - 1,
            "-",
            label=r"model $\beta_y$",
        )
        xlabel("s [m]")
        ylabel(r"$\Delta\beta$ [m]")

    def plot_dmux(m):
        muA = rad2deg(diff(unwrap(m.mux[m.xidxs] * 2 * pi)))
        muB = rad2deg(diff(unwrap(m.phase[m.xidxs])))
        spos = m.spos[m.xidxs][1:]
        res = m.res[m.xidxs]
        yerr = sqrt(res[1:] + res[:-1])
        pl.plot(spos, muA, "b-", label=r"model $\Delta\mu_x$")
        pl.plot(spos, muB, "r-", label=r"meas $\Delta\mu_x$")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def plot_dmuy(m):
        muA = rad2deg(diff(unwrap(m.muy[m.yidxs] * 2 * pi)))
        muB = rad2deg(diff(unwrap(m.phase[m.yidxs])))
        spos = m.spos[m.yidxs][1:]
        res = m.res[m.yidxs]
        pl.plot(spos, muA, "b-", label=r"model $\Delta\mu_y$")
        pl.plot(spos, muB, "r-", label=r"meas $\Delta\mu_y$")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def plot_dmux_beat(m):
        muA = rad2deg(diff(unwrap(m.mux[m.xidxs] * 2 * pi)))
        muB = rad2deg(diff(unwrap(m.phase[m.xidxs])))
        spos = m.spos[m.xidxs][1:]
        res = m.res[m.xidxs]
        yerr = sqrt(res[1:] + res[:-1])
        pl.plot(spos, muB - muA, "-", label=r"$\Delta\mu_x$ model-meas")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def plot_dmuy_beat(m):
        muA = rad2deg(diff(unwrap(m.muy[m.yidxs] * 2 * pi)))
        muB = rad2deg(diff(unwrap(m.phase[m.yidxs])))
        spos = m.spos[m.yidxs][1:]
        res = m.res[m.yidxs]
        pl.plot(spos, muB - muA, "-", label=r"$\Delta\mu_y$ model-meas")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def mmux(m):
        return cumsum((diff(unwrap(m.phase[m.xidxs]))) / 2 / pi)

    def mmuy(m):
        return cumsum((diff(unwrap(m.phase[m.yidxs]))) / 2 / pi)

    def plot_mux(m):
        delt = m.mux[m.xidxs][36] - m.mmux()[36]
        spos = m.spos[m.xidxs]
        pl.plot(spos, m.mux[m.xidxs], label=r"model $\mu_x$")
        mub = m.mmux()
        mub[36:] += delt
        pl.plot(spos[:-1], mub, label=r"meas $\mu_x$")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def plot_muy(m):
        delt = m.muy[m.yidxs][36] - m.mmuy()[36]
        spos = m.spos[m.yidxs]
        pl.plot(spos, m.muy[m.yidxs], label=r"model $\mu_y$")
        mub = m.mmuy()
        mub[36:] += delt
        pl.plot(spos[:-1], mub, label=r"meas $\mu_y$")
        xlabel("s [m]")
        ylabel(r"$\mu$ [deg]")
        legend()

    def plot_beta(m):
        m.plot_betx()
        m.plot_bety()
        legend()

    def plot_beta_beat(m):
        m.plot_betx_beat()
        m.plot_bety_beat()
        legend()


from . import rdmdate
from . import h5obj


class LHCBBQData(object):
    def __init__(self, filename):
        self.filename = filename
        self.__dict__.update(h5obj.load(self.filename, load_dataset=False))

    def find_idx(self, date):
        index = self.index[:]
        timestamp = index["timestamp"]
        idx = where(timestamp >= rdmdate.parsedate(date))[0][0] - 1
        ts = rdmdate.dumpdate(timestamp[idx])
        print("Found dataset %s" % ts)
        return index[idx]


from .harmonic_fit import *


def check_relrefmeas(vname, t1, t2):
    ldb = pytimber.LoggingDB()
    data = ldb.get(vname, t1, t2, scale="1 SECOND REPEAT", conf="mdb.conf")
    for i in range(0, len(vname), 2):
        print(vname[i])
        tm, vm = data[i]
        tr, vr = data[i + 1]
        vnorm = max(abs(vr))
        vnorm = 1
        pl.plot(tm - tm[0], (vr - vm) / vnorm, label=vname[i][:-7])
    pl.legend()
    return data


def check_refmeas(vname, t1, t2):
    ldb = pytimber.LoggingDB()
    data = ldb.get(vname, t1, t2, conf="mdb.conf")
    for i in range(0, len(vname), 2):
        print(vname[i])
        tm, vm = data[i]
        tr, vr = data[i + 1]
        pl.plot(tm - tm[0], vm, ".", label=vname[i])
        pl.plot(tr - tm[0], vr, "-", label=vname[i + 1])
    pl.legend()
    return data


def check_reflocal(vname, t1, t2, mask="data_ats_md3/power_conv2/%s.tsv.gz"):
    ldb = pytimber.LoggingDB()
    fnames = [mask % v for v in vname]
    data = ldb.get(fnames, t1=t1, t2=t2)
    for i in range(0, len(vname), 2):
        print(vname[i])
        tm, vm = data[i]
        tr, vr = data[i + 1]
        pl.plot(tm, vm, ".", label=vname[i])
        pl.plot(tr, vr, "-", label=vname[i + 1])
    pl.legend()
    return data


# data=check_refmeas(mq,'2011-06-29 04:30:00.000','2011-06-29 04:42:00.000')
#
#
# data=check_relrefmeas(mb,'2011-06-29 04:30:00.000','2011-06-29 04:42:00.000')
#
#
# figure()
# data=check_refmeas(mq,   '2011-06-29 06:00:00.000','2011-06-29 06:20:00.000')
# twinx()
# data=check_relrefmeas(mq,'2011-06-29 06:00:00.000','2011-06-29 06:20:00.000')
#
#
# data=check_relrefmeas(ir1,'2011-06-29 06:56:00.000','2011-06-29 07:20:00.000')
#
# data=check_refmeas(gen_both('.*RTQX.*L1.*'),'2011-06-29 06:56:00.000','2011-06-29 07:20:00.000')
