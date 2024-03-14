from table import stable
from plotfunc import plotfunc
import transport
from numpy import sqrt


class opticstable(stable, plotfunc):
    def maxbetx(f):
        return f.betx + f.alfx**2 / f.betx / abs(f.kn1l / f.l)

    def maxbety(f):
        return f.bety + f.alfy**2 / f.bety / abs(f.kn1l / f.l)

    def chromx(f):
        return -sum(f.kn1l * f.betx) / 4 / pi

    def chromy(f):
        return sum(f.kn1l * f.bety) / 4 / pi

    def ndx(t):
        return t.dx / sqrt(t.betx)

    def ndpx(t):
        return t.dpx * sqrt(t.betx) + t.dx / sqrt(t.betx) * t.alfx

    def alphac(t):
        return sum(t("dx*kn0l")) / sum(t.l)

    def gammatr(t):
        af = t.alphac()
        if af > 0:
            return sqrt(1 / af)
        else:
            return -sqrt(-1 / af)
