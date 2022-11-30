from numpy import *
from numpy.fft import rfft
from scipy.optimize import brent
from scipy.optimize import fmin_cg, fmin


def _fft_max(fv1):
    idx1 = abs(fv1).argmax()
    a11 = abs(fv1[idx1]) * 2
    p11 = angle(fv1[idx1])
    return idx1, a11, p11


def maxharm_fft(v1, fmin=0, fmax=0.5):
    n = float(len(v1))
    # fftmax 1
    fv1 = rfft(v1)
    fv1[: int(fmin * n)] = 0
    fv1[int(fmax * n) :] = 0
    idx1, a11, p11 = _fft_max(fv1)
    return idx1 / n, a11 / n, p11


def maxharm_coupled_fft(v1, v2):
    n1 = float(len(v1))
    n2 = float(len(v2))
    # fftmax 1
    fv1 = rfft(v1)
    co1 = fv1[0] / n1
    idx1, a11, p11 = _fft_max(fv1)
    ff1 = idx1 / n1
    a11 /= n1
    # fftmax 2
    fv2 = rfft(v2)
    co2 = fv2[0] / n2
    idx2, a22, p22 = _fft_max(fv2)
    ff2 = idx1 / n1
    a22 /= n2
    # cross terms
    a12 = abs(fv1[idx2]) * 2 / n1
    p12 = angle(fv1[idx2])
    a21 = abs(fv2[idx1]) * 2 / n2
    p21 = angle(fv2[idx1])
    return co1, co2, ff1, ff2, a11, a12, a21, a22, p11, p12, p21, p22


def _maxfreq_ftomin(ff, v, t):
    return -abs(sum(v * exp(2j * pi * ff * t)))


def _fit_maxharm_brent(v, fmax, t):
    ns = float(len(v))
    b = (fmax - 1.0 / ns, fmax, fmax + 1.0 / ns)
    try:
        fmax = brent(_maxfreq_ftomin, args=(v, t), brack=b, tol=1e-12)
    except ValueError:
        print("maxharm brent error")
        # print [(bb,_maxfreq_ftomin(bb,v,t)) for bb in b]
    coef = sum(v * exp(-2j * pi * fmax * t)) / ns
    a, p = abs(coef) * 2, angle(coef)
    res = sqrt(sum((v - a * cos(2 * pi * fmax * t + p)) ** 2) / sum(v**2)) / ns
    return fmax, a, p, res


def maxharm_brent(v, fmin=0, fmax=0.5):
    t = arange(len(v), dtype=float)
    ff, a, p = maxharm_fft(v, fmin=fmin, fmax=fmax)
    ff, a, p, res = _fit_maxharm_brent(v, ff, t)
    return ff, a, p, res


def _fit_maxharm_fmin(v, fmax, t):
    ns = float(len(v))
    fmax = fmin(_maxfreq_ftomin, fmax, args=(v, t), disp=False, xtol=1e-10)
    coef = sum(v * exp(-2j * pi * fmax * t)) / ns
    return fmax, abs(coef) * 2, angle(coef)


def maxharm_fftmax(v):
    t = arange(len(v), dtype=float)
    ff, a, p = maxharm_fft(v)
    return _fit_maxharm_fmin(v, ff, t)


def maxharm_lsq(v1):
    ns = float(len(v1))
    ff1, a11, p11 = maxharm_fft(v1)
    x = ff1, a11, p11
    t = arange(ns)

    def ftomin(x):
        ff1, a11, p11 = x
        vv1 = a11 * cos(2 * pi * ff1 * t + p11)
        return sum((vv1 - v1) ** 2)

    ff1, a11, p11 = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin((ff1, a11, p11))
    return ff1, a11, p11, res


def maxharm_lsq2(v1, v2):
    ns = float(len(v1))
    ff1, a11, p11 = maxharm_fft(v1)
    ff2, a21, p21 = maxharm_fft(v2)
    x = ff1, a11, p11, a21, p21
    t = arange(ns)

    def ftomin(x):
        ff1, a11, p11, a21, p21 = x
        vv1 = a11 * cos(2 * pi * ff1 * t + p11)
        vv2 = a21 * cos(2 * pi * ff1 * t + p21)
        return sum((vv1 - v1) ** 2 + (vv2 - v2) ** 2)

    ff1, a11, p11, a21, p21 = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin((ff1, a11, p11, a21, p21))
    return ff1, a11, p11, a21, p21, res


unit = lambda n: 1.0
sinw = lambda n: sin(linspace(0, pi, n))
sin2w = lambda n: sin(linspace(0, pi, n)) ** 2
sin4w = lambda n: sin(linspace(0, pi, n)) ** 4
sin8w = lambda n: sin(linspace(0, pi, n)) ** 8


def _get_harm_coef(t, v, fmax):
    coef = sum(v * exp(-2j * pi * fmax * t)) / len(t)
    return abs(coef) * 2, angle(coef)


def get_harm_coef(v, f):
    ns = float(len(v))
    t = arange(ns)
    return _get_harm_coef(t, v, f)


def fit_coupled_lsq(v1, v2):
    ns = float(len(v1))
    co1 = mean(v1)
    co2 = mean(v2)
    ff1, a11, p11, res = maxharm_lsq(v1 - co1)
    ff2, a22, p22, res = maxharm_lsq(v2 - co2)
    t = arange(ns)
    a12, p12 = _get_harm_coef(t, v1, ff2)
    a21, p21 = _get_harm_coef(t, v2, ff1)
    x = a12, p12

    def ftomin(x):
        a12, p12 = x
        vv1 = mk_harmreal_vec(t, [0, ff1, ff2], [co1, a11, a12], [0, p11, p12])
        res = sum((vv1 - v1) ** 2)
        return res

    x = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin(x)
    a21, p21 = x
    x = a21, p21

    def ftomin(x):
        a21, p21 = x
        vv2 = mk_harmreal_vec(t, [0, ff1, ff2], [co1, a21, a22], [0, p21, p22])
        res = sum((vv2 - v2) ** 2)
        return res

    x = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    a21, p21 = x
    res += ftomin(x)
    return co1, co2, ff1, ff2, a11, a12, a21, a22, p11, p12, p21, p22, res


def fit_coupled_lsq2(v1, v2):
    ns = float(len(v1))
    co1 = mean(v1)
    co2 = mean(v2)
    ff1, a11, p11, res = maxharm_lsq(v1 - co1)
    ff2, a22, p22, res = maxharm_lsq(v2 - co2)
    t = arange(ns)
    a12, p12 = _get_harm_coef(t, v1, ff2)
    a21, p21 = _get_harm_coef(t, v2, ff1)
    x = co1, co2, ff1, ff2, a11, a12, a21, a22, p11, p12, p21, p22

    def ftomin(x):
        co1, co2, ff1, ff2, a11, a12, a21, a22, p11, p12, p21, p22 = x
        vv1 = (
            co1 + a11 * cos(2 * pi * ff1 * t + p11) + a12 * cos(2 * pi * ff2 * t + p12)
        )
        vv2 = (
            co2 + a21 * cos(2 * pi * ff1 * t + p21) + a22 * cos(2 * pi * ff2 * t + p22)
        )
        res = sum((vv1 - v1) ** 2) / ns / sum(v1**2) + sum(
            (vv2 - v2) ** 2
        ) / ns / sum(v2**2)
        return res

    x = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin(x)
    return co1, co2, ff1, ff2, a11, a12, a21, a22, p11, p12, p21, p22, res


def fit_single_lsq(v1):
    ns = float(len(v1))
    co1 = mean(v1)
    ff1, a11, p11, res = maxharm_lsq(v1 - co1)
    return co1, ff1, a11, p11, res


def mk_harmreal_vec(t, freq, ampl, phase):
    s = t * 0
    for ff, a, p in zip(freq, ampl, phase):
        s += a * cos(2 * pi * ff * t + p)
    return s


def harm_comp_rel(s, n=100, eps=1e-5, refit=False):
    s = s.copy()
    ns = float(len(s))
    t = arange(ns)
    c = mean(s)
    out = [[0, c, 0]]
    s -= c
    for i in range(n):
        error = sum(abs(s) ** 2)
        if error < eps:
            break
        ff, a, p = maxharm_fft(s)
        ff, a, p = _fit_maxharm_brent(s, ff, t)
        if refit is True:
            ff, a, p, res = maxharm_lsq(ff, a, p)
        sfit = c * cos(2 * pi * ff * t + p)
        out.append([f, a, p])
        s -= sfit
    return list(zip(*out))


def mk_coupled_vec(t, co1, co2, f1, f2, a11, a12, a21, a22, p11, p12, p21, p22):
    t1 = 2 * pi * f1 * t
    t2 = 2 * pi * f2 * t
    vv1 = co1 + a11 * cos(t1 + p11) + a12 * cos(t2 + p12)
    vv2 = co2 + a21 * cos(t1 + p21) + a22 * cos(t2 + p22)
    return vv1, vv2


if __name__ == "__main__":
    from numpy import random

    t = arange(3000, dtype=float)
    noise1 = 0.1 * (random.rand(3000) - 0.5)
    noise2 = 0.1 * (random.rand(3000) - 0.5)
    v1 = mk_harmreal_vec(t, [0, 0.31, 0.28], [2, 1, 0.1], [0, 0.3, 0.4]) + noise1
    v2 = mk_harmreal_vec(t, [0, 0.31, 0.28], [3, 0.1, 1], [0, 0.2, 0.1]) + noise2
    print(maxharm_fft(v1 - mean(v1)))
    print(maxharm_brent(v1 - mean(v1)))
    print(maxharm_lsq(v1 - mean(v1)))
    print(fit_single_lsq(v1))
    print(fit_coupled_lsq(v1, v2))


def maxharm_lsq(v1):
    ns = float(len(v1))
    ff1, a11, p11 = maxharm_fft(v1)
    x = ff1, a11, p11
    t = arange(ns)

    def ftomin(x):
        ff1, a11, p11 = x
        vv1 = a11 * cos(2 * pi * ff1 * t + p11)
        return sum((vv1 - v1) ** 2)

    ff1, a11, p11 = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin((ff1, a11, p11))
    return ff1, a11, p11, res


def f_decay_lin(t, ff, a0, a1, p, m0):
    return m0 + (a0 * exp(t / a1)) * cos(2 * pi * ff * t + p)


def fit_decay(v, fa, fb):
    ns = float(len(v))
    ff, a0, p = maxharm_fft(v, fa, fb)
    m0 = v.mean()
    x = ff, a0, -1.0, p, m0
    t = arange(ns)

    def ftomin(x):
        ff, a0, a1, p, m0 = x
        vv = f_decay_lin(t, ff, a0, a1, p, m0)
        res = sum((vv - v) ** 2) / ns / sum(v) ** 2
        return res

    ff, a0, a1, p, m0 = fmin(ftomin, x, xtol=1e-10, ftol=1e-10, disp=False)
    res = ftomin((ff, a0, a1, p, m0))
    return ff, a0, a1, p, m0, res
