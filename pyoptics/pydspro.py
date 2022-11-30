from matplotlib.pyplot import subplot, ylabel, xlabel, plot, grid, imshow, title
import matplotlib.pyplot as _p
from numpy import *
from scipy.special import erf
from numpy.fft import *

# from exprtotex import e2t


def c2db(c):
    """complex signal to decibel"""
    return 20.0 * log10(abs(c))


def db2c(d):
    """decibel to complex"""
    return 10 ** (d / 20.0)


def a2c(a):
    """angle to complex"""
    return exp(1j * pi * a / 180.0)


def c2a(c):
    """complex to angle"""
    return angle(c) / pi * 180.0


def c2ua(c):
    return unwrap(c2a(c), discont=180.0)


def gdelay(f, c):
    df = f[1] - f[0]
    return r_[0.0, diff(unwrap(angle(c))) / 2 / pi / df]


def resample(y, newsize):
    oldsize = y.size
    tempsize = oldsize * newsize
    tempy = zeros(tempsize, dtype=y.dtype)
    newy = zeros(newsize, dtype=y.dtype)
    for i in range(newsize):
        tempy[i::newsize] = y
    for i in range(oldsize):
        newy += tempy[i::oldsize]
    newy /= oldsize
    return newy


c = 299792458.0  # m/s
epsilon0 = 8.8541878176e-12  # F/m
mu0 = pi * 4e-7  # H/m
z0 = 376.730313462  # Ohm


class plotdelay:
    def __init__(self, f, v):
        self.d1, self.d0 = 0.0, 0.0
        self.f = f
        self.v = v
        (self.p,) = plot(f, c2a(v), label="%9.4e" % self.d1)
        legend()

    def redraw(self):
        self.p.set_label("%9.4e" % self.d1)
        self.p.set_ydata(c2a(self.v * delay(self.f, [self.d0, self.d1])))
        legend()

    def __mul__(self, c):
        self.d1 *= c
        self.redraw()

    def __div__(self, c):
        self.d1 /= c
        self.redraw()

    def set(self, c):
        self.d1 = c
        self.redraw()


def plotdba(f, v, title="", label="", wrap=True):
    """plot db and angle"""
    p1 = subplot(211)
    plot(f, c2db(v))
    if title:
        _p.title(title)
    ylabel("abs [dB]")
    grid(True)
    p2 = subplot(212, sharex=p1)
    if wrap:
        plot(f, c2a(v), label=label)
    else:
        plot(f, c2ua(v), label=label)
    ylabel("angle [degree]")
    xlabel("freq [Hz]")
    grid(True)
    return p1, p2


def plotpolar(f, v, wrap=True):
    """plot db and angle"""
    p1 = subplot(211)
    plot(f, v)
    ylabel("abs")
    grid(True)
    p2 = subplot(212)
    a = angle(v)
    if not wrap:
        a = unwrap(a, discont=pi / 2)
    plot(f, a)
    ylabel("angle")
    xlabel("freq [Hz]")
    grid(True)
    return p1, p2


def gauss(t, mu, sigma):
    return 1.0 / sigma / sqrt(2 * pi) * exp(-0.5 * ((t - mu) / sigma) ** 2)


def rect(t, mu, sigma):
    return 1.0 * (abs(t - mu) < sigma)


def phi(t):
    return 0.5 * (1 + erf(t / sqrt(2)))


def delay(f, dt):
    """Return a delay in power of frequency
    2j*pi*dt[0]+2j*pi*dt[1]*f+...
    """
    if not hasattr(dt, "__len__"):
        dt = [dt]
    a = zeros(len(f))
    nf = 1
    for i in dt:
        a += i * nf
        nf *= f
    return exp(2j * pi * a)


def sumsq(x):
    return sum(abs(x) ** 2)


def air_cable(f, a1, a2):
    """air cable"""
    return exp(-a1 * sqrt(2j * f) - a2 * f)


def _time(n, fs):
    fs = float(fs)
    return arange(0.0, n) / fs


def _freq(n, fs):
    fs = float(fs)
    #  return arange(0,fs/2+fs/n,fs/n)
    return arange(0.0, n / 2 + 1) * fs / n


def decay(t, st, tau):
    s = exp(-(t - st) / tau)
    s[t < st] = 0
    return s


def t2f(t):
    n = len(t) / 2 + 1
    fs = 1 / (t[1] - t[0])
    return arange(0.0, n) * (fs / len(t))


def f2t(f):
    n = (len(f) - 1) * 2
    fs = 2 * f[-1]
    return arange(0.0, n) / fs


def ff(n):
    return arange(0, 1, 1.0 / n)


def rff(n):
    return arange(0, 0.5 + 0.5 / n, 1.0 / n)


def tt(n):
    return arange(n)


def nextp2(n):
    return int(2 ** ceil(log2(n)))


def padr(f, n=None):
    if n is None:
        n = nextp2(len(f))
    return r_[f, zeros(n - len(f))]


def padf(f, n=None):
    if n is None:
        n = nextp2(len(f))
    n = n * 2
    fs = n * (f[1] - f[0])
    return freq(n, fs)


def padc(f, n=None):
    if n is None:
        n = nextp2(len(f))
    nz = n - len(f)
    if nz % 2 == 0:
        n1 = nz // 2
        n2 = n1
    else:
        n1 = nz // 2
        n2 = n1 + 1
    return r_[zeros(n1), f, zeros(n2)]


def pad_right(a, l):
    return r_[a, zeros(l - len(a), dtype=a.dtype)]


def pad_left(a, l):
    return r_[zeros(l - len(a), dtype=a.dtype), a]


def pad_center(a, l):
    ll = (l - len(a)) / 2
    lr = l - len(a) - ll
    return r_[zeros(ll, dtype=a.dtype), a, zeros(lr, dtype=a.dtype)]


def pad(t, s):
    ss = zeros(len(t), dtype=s.dtype)
    start = (len(t) - len(s)) / 2
    ss[start : start + len(s)] = s
    return ss


def rect(f, a):
    s = ones(len(f), dtype=float)
    w = f > a
    s[w] = 0
    return s


def hann(f, a):
    r = a / f[-1]
    r = r % 1
    ns = floor((len(f) + len(f) % 2 - 1) * r)
    s = zeros_like(f)
    #  print ns
    s[:ns] = hanning(2 * ns)[ns:]
    return s


def hann2(f, l, h):
    s = ones_like(f)
    il = f < l
    l = sum(il)
    s[:l] = hanning(2 * l)[:l]
    il = f < h
    l = sum(il)
    s[-l:] = hanning(2 * l)[l:]
    return s


def c2gd(a):
    a = r_[0, diff(c2a(a).copy())]
    a[a > 350] -= 360
    return a


def upmethod(cls):
    src = []
    s = None
    while s != "":
        s = input()
        src.append(s[2:])
    src = "\n".join(src)
    c = compile(src, "exec", "exec")
    eval(c)
    name = c.co_names[0]
    setattr(cls, name, locals()[name])


def expcouplertf(f, l=0.40, a=log(0.03 / 0.002), v=c):
    t = f2t(f)
    s = 0 * t
    tt = t < 2 * l / v
    s[tt] = 0.5 * exp(-v * t[tt] * a / (2 * l))
    s = s / sum(s)
    tf = 2j * pi * f * rfft(s)
    tf[0] = 1e12
    return tf


def window(start, stop, l, f):
    return r_[zeros(start), f(stop - start), zeros(l - stop)]


def plotimg(f, s, **args):
    imshow(s, origin="bottom", aspect="auto", extent=[0, f[-1], 0, s.shape[0]], **args)


def edgedetect(v, ma, tr):
    w = convolve(abs(v), 1.0 / ma * ones(ma), mode="same")
    ww = w.copy()
    w[ww > tr] = 1
    w[ww <= tr] = 0
    #  w=r_[w[ma/2:],w[:ma/2]]
    w = diff(w)
    return w


def edgeget(v, w, wl, bn):
    start = where(w == 1)[0]
    stop = where(w == -1)[0]
    w = c_[start, stop]
    c = (w[bn, 0] + w[bn, 1]) / 2
    return v[:, c - wl / 2 : c + wl / 2]


def movavg2(x, n=100):
    """moving average with rect"""
    return convolve(ones(n, dtype=float) / n, x, mode="same")


def movavg3(x, n=100, reverse=False):
    """moving average with rect, fixside"""
    x = convolve(ones(n, dtype=float) / n, x, mode="same")
    x[: n / 2] = x[n / 2]
    return x


def rng(x, a, b):
    "return (x<b) & (x>a)"
    return (x < b) & (x > a)


def dBm2Vpp(x, z=50.0):
    return sqrt(8 * z * 10 ** (x * 0.1 - 3))


def Vpp2dBm(x, z=50.0):
    return 10 * log10(x**2 / z * 125)


def h5print(x, pk=""):
    print(pk, "->", end=" ")
    if hasattr(x, "shape") and hasattr(x, "dtype"):
        print("%s[%s]" % (x.dtype.str, ",".join(map(str, x.shape))), end=" ")
    elif type(x) in [int, float, int]:
        print(x, end=" ")
    elif type(x) in [str]:
        if len(x) > 20:
            print("%s..%s" % (x[:8], x[-8:]), end=" ")
        else:
            print(x)
    else:
        print(type(x).__name__, end=" ")
    print()
    if hasattr(x, "keys"):
        if "h5py" in str(x.__class__):
            for k in list(x.keys()):
                if len(pk) == 0:
                    npk = "['%s']" % k
                else:
                    npk = pk[:-2] + "/%s']" % k
                h5print(x[k], pk=npk)
        else:
            for k in list(x.keys()):
                npk = pk + "['%s']" % k
                h5print(x[k], pk=npk)
    if hasattr(x, "attrs"):
        for k in list(x.attrs.keys()):
            npk = pk + ".attrs['%s']" % k
            try:
                h5print(x.attrs[k], pk=npk)
            except TypeError:
                print("TypeError")


def reduceAx(v, n, axis=0):
    sl = [slice(None)] * len(v.shape)
    sl[axis] = slice(0, None, n)
    vnew = v[tuple(sl)].copy()
    for i in range(1, n):
        sl[axis] = slice(i, None, n)
        vnew += v[tuple(sl)]
    vnew /= n
    return vnew


def smoothListGaussian(list, degree=5):
    window = degree * 2 - 1
    padl = window / 2
    list = r_[zeros(padl), list, zeros(window - padl)]
    weight = array([1.0] * window)
    weightGauss = []
    for i in range(window):
        i = i - degree + 1
        frac = i / float(window)
        gauss = 1 / (exp((4 * (frac)) ** 2))
        weightGauss.append(gauss)
    weight = array(weightGauss) * weight
    smoothed = [0.0] * (len(list) - window)
    for i in range(len(smoothed)):
        smoothed[i] = sum(array(list[i : i + window]) * weight) / sum(weight)
    return smoothed


def peak_detect(y, wdw=10):
    """Return an array of local maxima within a window"""
    leny = len(y)
    newy = r_[zeros(wdw), y, zeros(wdw)]
    vl = [roll(newy, i) for i in range(-wdw, wdw)]
    yy = vstack(vl).max(axis=0)[wdw:-wdw]
    return yy
