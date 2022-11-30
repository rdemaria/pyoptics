import os

import matplotlib.pyplot as pl
import matplotlib.delaunay
from convexhull import convex_hull

from numpy import *
import tfsdata
from pyoptics import optics


def find_res_xcross(m, n, q, xs, y1, y2, out):
    #  print 'x=%d'% (xs),
    if n != 0:
        m, n, q, xs, y1, y2 = map(float, (m, n, q, xs, y1, y2))
        ys = (q - m * xs) / n
        #    print 'ys=%g'%ys,
        #    print '%g<=ys<=%g'%(y1,y2),
        if ys >= y1 and ys <= y2:
            #      print [xs,ys],
            out.append((xs, ys))


#  print


def find_res_ycross(m, n, q, ys, x1, x2, out):
    #  print 'y=%d'% (ys),
    if m != 0:
        m, n, q, ys, y1, y2 = map(float, (m, n, q, ys, x1, x2))
        xs = (q - n * ys) / m
        #    print 'xs=%g'%xs,
        #    print '%g<=xs<=%g'%(x1,x2),
        if xs >= x1 and xs <= x2:
            #      print [xs,ys]
            out.append((xs, ys))


#  print


def get_res_box(m, n, a=0, b=1, c=0, d=1):
    """m,n,q resonance integers
    a,b,c,d box parameters
    """
    order = int(_n.ceil(abs(m) * max(abs(a), abs(b)) + abs(n) * max(abs(c), abs(d))))
    out = []
    for q in range(-order, +order + 1):
        points = []
        #    print '%2d*x+%2d*y=%2d' % (m,n,q)
        find_res_xcross(m, n, q, a, c, d, points)
        find_res_xcross(m, n, q, b, c, d, points)
        find_res_ycross(m, n, q, c, a, b, points)
        find_res_ycross(m, n, q, d, a, b, points)
        points = list(set(points))
        if len(points) > 1:
            out.append(points)
    return out


def getmn(o, s="b"):
    out = []
    for m in range(0, o + 1):
        n = o - m
        if s == "b" and n % 2 == 0:
            out.append((m, n))
            if n > 0:
                out.append((m, -n))
        elif s == "a" and n % 2 == 1:
            out.append((m, n))
            if n > 0:
                out.append((m, -n))
    return out


def plot_res_box(m, n, a=0, b=1, c=0, d=1, color="k"):
    points = get_res_box(m, n, a, b, c, d)
    for c in points:
        x, y = zip(*c)
        pl.plot(x, y, "%s" % color)


def plot_res_order_box(o, a=0, b=1, c=0, d=1, c1="k", c2="r"):
    for m, n in getmn(o, "b"):
        # print 'b%s: m=%d n=%d'%(o,m,n)
        plot_res_box(m, n, a=a, b=b, c=c, d=d, color=c1)
    for m, n in getmn(o, "a"):
        # print 'a%s: m=%d n=%d'%(o,m,n)
        plot_res_box(m, n, a=a, b=b, c=c, d=d, color=c2)


def plot_res(m, n, color="k"):
    a, b = pl.xlim()
    c, d = pl.ylim()
    points = get_res_box(m, n, a, b, c, d)
    for c in points:
        x, y = zip(*c)
        pl.plot(x, y, "-%s" % color)


def plot_res_order(o, c1="k", c2="r"):
    a, b = pl.xlim()
    c, d = pl.ylim()
    for m, n in getmn(o, "b"):
        # print 'b%s: m=%d n=%d'%(o,m,n)
        plot_res_box(m, n, a=a, b=b, c=c, d=d, color=c1)
    for m, n in getmn(o, "a"):
        # print 'a%s: m=%d n=%d'%(o,m,n)
        plot_res_box(m, n, a=a, b=b, c=c, d=d, color=c2)


def mkranges(nsigma=12, nangles=7):
    ranges = []
    for i in range(nangles):
        ranges.append([0, i])
    for i in range(nangles):
        ranges.append(slice(1 + i, nangles * nsigma + 1, nangles))
    for i in range(nsigma):
        ranges.append(slice(1 + nangles * i, 1 + nangles * (i + 1)))
    return ranges


mycolors = list("rcgmb")


def colorrotate():
    c = mycolors.pop(0)
    mycolors.append(c)
    return c


def change_alpha(ll, alpha):
    for a in ll:
        hasasttr(a, "set_alpha") and a.set_alpha(alpha) or change_alpha(a, alpha)


def Dq2(b2, Bx, Jx, By, Jy):
    #  b2=k1l
    Dqx = b2 * Bx / (4.0 * pi)
    Dqy = -b2 * By / (4.0 * pi)
    return Dqx, Dqy


def Dq4(b4, Bx, Jx, By, Jy):
    #  b4=k3l/6.
    Dqx = b4 * 3 * Bx * (Bx * Jx - 2 * By * Jy) / (8.0 * pi)
    Dqy = b4 * 3 * By * (By * Jy - 2 * Bx * Jx) / (8.0 * pi)
    return Dqx, Dqy


def Dq6(b6, Bx, Jx, By, Jy):
    #  b6=k5l/120.
    Dqx = (
        b6
        * 5
        * Bx
        * (Bx**2 * Jx**2 - 6 * Bx * By * Jx * Jy + 3 * By**2 * Jy**2)
        / (8.0 * pi)
    )
    Dqy = (
        -b6
        * 5
        * By
        * (3 * Bx**2 * Jx**2 - 6 * Bx * By * Jx * Jy + By**2 * Jy**2)
        / (8.0 * pi)
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
        / (32.0 * pi)
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
        / (32.0 * pi)
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
        / (32.0 * pi)
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
        / (32.0 * pi)
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
        / (64.0 * pi)
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
        / (64.0 * pi)
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
        / (64.0 * pi)
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
        / (64.0 * Pi)
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
        / (512.0 * pi)
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
        / (512.0 * pi)
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


def kTob(kn, ks, x, y):
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


class Footprint(object):
    def plot_grid(t, nsigma=None, lw=1):
        if nsigma is None:
            nsigma = t.nsigma
        nangles = t.nangles
        ranges = mkranges(nsigma, nangles)
        for i in ranges:
            if hasattr(i, "step"):
                lw = i.step == nangles and i.start / 2.0 or 1
            pl.plot(t.x[i], t.y[i], "-k", lw=lw)

    def plot_footprint(
        t, nsigma=12, nangles=7, wp=(0.28, 0.31), spread=0.01, label=None, color=None
    ):
        ranges = mkranges(nsigma, nangles)
        lw = 1
        out = []
        lbl = True
        if label is None:
            label = t.label
        if color is None:
            color = colorrotate()
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
            plot(t.tunx[i], t.tuny[i])

    def reshape(self):
        """return tunes in [sigma,angles]"""
        qx = self.tunx[1:].reshape(self.nsigma, self.nangles)
        qy = self.tuny[1:].reshape(self.nsigma, self.nangles)
        return qx, qy


class FootTrack(Footprint):
    def __init__(self, base, nangles=7, nsigma=12):
        self.base = base
        self.label = base.replace("_", " ")
        dynapfn = os.path.join(base, "dynaptune.tfs")
        t = tfsdata.open(dynapfn)
        self.tunx = t["tunx"]
        self.tuny = t["tuny"]
        self.tx = t["x"]
        self.ty = t["y"]
        self.nangles = nangles
        self.nsigma = nsigma

    def mkComp(self):
        return FootComp(self.base, nangles=self.nangles)


class FootComp(Footprint):
    def __init__(self, base, nsigma=12, nangles=7, wp=(0.28, 0.31), mcx=False):
        self.base = base
        self.label = base.replace("_", " ") + " comp"
        self.tw = optics.open(os.path.join(base, "optics0_inser.mad"))
        self.et = optics.open(os.path.join(base, "tripD1D2.errors"))
        if mcx:
            self.mcx = os.path.join(base, "MCX_setting.mad")
        else:
            self.mcx = False
        self.nsigma = nsigma
        self.nangles = nangles
        self.x, self.y = self.mk_grid()
        self.tunx = self.x * 0 + wp[0]
        self.tuny = self.y * 0 + wp[1]
        self.ex = self.tw.header["ex"]
        self.ey = self.tw.header["ex"]
        self.mk_detuning()

    def change_emit(self, emit_n=3.75e-6, energy=450, pmass=0.938):
        self.ex = emit_n / energy * pmass
        self.ey = emit_n / energy * pmass

    def mk_grid(self):
        small = 0.05
        big = sqrt(1.0 - small**2)
        n = 1
        m = 0
        # sigma angle multiplier
        x = [small]
        y = [small]
        while n <= self.nsigma:
            angle = 90.0 / (self.nangles - 1) * m * pi / 180
            if m == 0:
                xs = n * big
                ys = n * small
            elif m == self.nangles - 1:
                xs = n * small
                ys = n * big
            else:
                xs = n * cos(angle)
                ys = n * sin(angle)
            m = m + 1
            if m == self.nangles:
                m = 0
                n = n + 1
            x.append(xs)
            y.append(ys)
        return array(x), array(y)

    def mk_detuning(self, mpole="k5l", feeddown=False):
        dqxl = {}
        dqyl = {}
        data = self.get_twiss_err(mpole)
        if self.mcx:
            if mpole == "k5l":
                data += self.get_mcx_setting("KCTX")
        dqxt = []
        dqyt = []
        i = 0
        for xs, ys in zip(self.x, self.y):
            dqxll = []
            dqyll = []
            for name, k5l, Bx, By, x, y in data:
                b6 = k5l / 120
                b4 = FeedDown([0, 0, 0, 0, 0, b6], [0, 0, 0, 0, 0, 0], x, y, 4)[0]
                # b2=FeedDown([0,0,0,0,0,b6],[0,0,0,0,0,0],x,y,2)[0]
                Jx = self.ex / 2 * xs**2
                Jy = self.ey / 2 * ys**2
                dqx6, dqy6 = Dq6(b6, Bx, Jx, By, Jy)
                dqx4, dqy4 = Dq4(b4, Bx, Jx, By, Jy)
                # dqx2,dqy2=Dq2(b2,Bx,Jx,By,Jy)
                dqxll.append(dqx4 + dqx6)
                dqyll.append(dqy4 + dqy6)
            self.tunx[i] += sum(dqxll)
            self.tuny[i] += sum(dqyll)
            i += 1
        return self

    def get_twiss_err(self, mpole):
        """err =e1.k5l"""
        e1 = self.et
        t1 = self.tw
        err = e1[mpole]
        idxe1 = where(abs(err) > 0)[0]
        names = e1.name[idxe1]
        out = []
        for name, k5l in zip(e1.name[idxe1], err[idxe1]):
            idxt1 = where(t1.name == name)[0][0]
            Bx, By = t1.betx[idxt1], t1.bety[idxt1]
            x, y = t1.x[idxt1], t1.y[idxt1]
            # Dx,Dy=t1.dx[idxt1],t1.dy[idxt1]
            out.append([name, k5l, Bx, By, x, y])
        return out

    def get_mcx_setting(self, fam):
        data = {}
        for line in open(self.mcx):
            if ":=" in line:
                ll = line.split()
                name = ll[0]
                if fam in name:
                    val = float(ll[2])
                    n = name
                    nname = "M" + n[1:4] + "." + n[-4] + n[-2:]
                    data[nname] = val
        t1 = self.tw
        out = []
        for name, k5l in data.items():
            if name in t1.name:
                idxt1 = where(t1.name == name)[0][0]
                Bx, By = t1.betx[idxt1], t1.bety[idxt1]
                x, y = t1.x[idxt1], t1.y[idxt1]
                out.append([name, k5l, Bx, By, x, y])
        return out


# f1=FootComp('temp')
# f2=FootTrack('temp')
