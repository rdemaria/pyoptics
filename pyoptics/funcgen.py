import numpy as np
import scipy.optimize


class PPLP(object):
    def __init__(self, m0, a1, t1, a2, t2, a3, t3, a4):
        self.m0 = m0
        self.a1 = a1
        self.t1 = t1
        self.a2 = a2
        self.t2 = t2
        self.a3 = a3
        self.t3 = t3
        self.a4 = a4

    def __repr__(self):
        m0 = self.m0
        a1 = self.a1
        t1 = self.t1
        a2 = self.a2
        t2 = self.t2
        a3 = self.a3
        t3 = self.t3
        a4 = self.a4
        return "PPLP(%s)" % (",".join(map(str, (m0, a1, t1, a2, t2, a3, t3, a4))))

    def fit_rate_energy(self, energy=7000e9, rate=5.7967e9):
        def ftosolve(x):
            self.t2 = x[0]
            self.t3 = x[1]
            tf = self.duration()
            return self.derivative(150) - rate, self(tf) - energy

        self.t2, self.t3 = scipy.optimize.fsolve(ftosolve, (self.t2, self.t3))
        return self

    def duration(self):
        m0 = self.m0
        a1 = self.a1
        t1 = self.t1
        a2 = self.a2
        t2 = self.t2
        a3 = self.a3
        t3 = self.t3
        a4 = self.a4
        d1 = a1 * t1
        d2 = d1 + a2 * (t2 - t1)
        d3 = d2 + a3 * (t3 - t2)
        return t3 - d3 / a4

    def __call__(self, tt):
        t0 = 0
        m0 = self.m0
        a1 = self.a1
        t1 = self.t1
        a2 = self.a2
        t2 = self.t2
        a3 = self.a3
        t3 = self.t3
        a4 = self.a4
        c0 = m0
        d0 = 0
        c1 = c0 + d0 * t1 + a1 / 2 * t1**2
        d1 = a1 * t1
        d2 = d1 + a2 * (t2 - t1)
        d3 = d2 + a3 * (t3 - t2)
        c2 = c1 + d1 * (t2 - t1) + a2 / 2 * (t2 - t1) ** 2
        c3 = c2 + d2 * (t3 - t2) + a3 / 2 * (t3 - t2) ** 2
        tf = t3 - d3 / a4
        c4 = c3 + d3 * (tf - t3) + a4 / 2 * (tf - t3) ** 2
        tt1 = tt - t1
        tt2 = tt - t2
        tt3 = tt - t3
        mom = (
            c0
            if (tt < t0)
            else c0 + a1 * tt**2 / 2
            if (t0 <= tt) and (tt < t1)
            else c1 + d1 * tt1 + a2 * tt1**2 / 2
            if (t1 <= tt) and (tt < t2)
            else c2 + d2 * tt2 + a3 * tt2**2 / 2
            if (t2 <= tt) and (tt < t3)
            else c3 + d3 * tt3 + a4 * tt3**2 / 2
            if (t3 <= tt) and (tt < tf)
            else c4
        )
        return mom

    def derivative(self, tt):
        m0 = self.m0
        a1 = self.a1
        t1 = self.t1
        a2 = self.a2
        t2 = self.t2
        a3 = self.a3
        t3 = self.t3
        a4 = self.a4
        c0 = m0
        d0 = 0
        c1 = c0 + d0 * t1 + a1 / 2 * t1**2
        d1 = a1 * t1
        d2 = d1 + a2 * (t2 - t1)
        d3 = d2 + a3 * (t3 - t2)
        tt1 = tt - t1
        tt2 = tt - t2
        tt3 = tt - t3
        tf = t3 - d3 / a4
        momp = (
            a1 * tt
            if (tt < t1)
            else d1 + a2 * tt1
            if (t1 <= tt) and (tt < t2)
            else d2 + a3 * tt2
            if (t2 <= tt) and (tt < t3)
            else d3 + a4 * tt3
            if (t3 <= tt) and (tt < tf)
            else 0
        )
        return momp


class PELP(object):
    def __init__(self, Ii, If, A, D, R, Te, Ti):
        self.Ii = Ii
        self.If = If
        self.A = A
        self.D = D
        self.R = R
        self.Te = Te
        self.Ti = Ti

    def duration(self):
        Ii = self.Ii
        If = self.If
        A = self.A
        D = self.D
        R = self.R
        Te = self.Te
        Ti = self.Ti
        iaTe = A / 2 * (Te - Ti) ** 2 + Ii
        iapTe = A * (Te - Ti)
        b = iapTe / iaTe
        a = iaTe * np.exp(-b * Te)
        tl = np.log(R / (a * b)) / b if Te != 0 else R / A + Ti
        il = a * np.exp(b * tl) if Te != 0 else A / 2 * (tl - Ti) ** 2 + Ii
        td = (If - il) / R + tl - R / (2 * D)
        tf = td + R / D
        return tf

    def __call__(self, time):
        Ii = self.Ii
        If = self.If
        A = self.A
        D = self.D
        R = self.R
        Te = self.Te
        Ti = self.Ti
        iaTe = A / 2 * (Te - Ti) ** 2 + Ii
        iapTe = A * (Te - Ti)
        b = iapTe / iaTe
        a = iaTe * np.exp(-b * Te)
        tl = np.log(R / (a * b)) / b if Te != 0 else R / A + Ti
        il = a * np.exp(b * tl) if Te != 0 else A / 2 * (tl - Ti) ** 2 + Ii
        td = (If - il) / R + tl - R / (2 * D)
        tf = td + R / D
        I = (
            A / 2 * (time - Ti) ** 2 + Ii
            if (Ti <= time) and (time < Te)
            else a * np.exp(b * time)
            if (Te <= time) and (time < tl)
            else R * (time - tl) + il
            if (tl <= time) and (time < td)
            else -D / 2 * (tf - time) ** 2 + If
            if (td <= time) and (time < tf)
            else If
        )
        return I


pelp7000 = PELP(450e9, 7000e9, 0.009, 0.02, 10.0, 325, 0)
pelp6500 = PELP(450e9, 6500e9, 0.009, 0.02, 10.0, 325, 0)
pplp6500 = PPLP(450e9, 0.0133e9, 38, 0.3e9, 55.6376, 0, 1078.429, -0.3e9)
pplp7000 = PPLP(450e9, 0.0133e9, 38, 0.3e9, 55.6376, 0, 1164.685, -0.3e9)
pplp6500.fit_rate_energy(6500e9)
pplp7000.fit_rate_energy(7000e9)
