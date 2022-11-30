from __future__ import division

from numpy import *
import matplotlib.pyplot as plt

try:
    from objdebug import ObjDebug
except ImportError:
    ObjDebug = object


class Quantity(object):
    def __init__(self, default, desc, fmt="%g"):
        self.desc = desc
        self.fmt = fmt
        self.default = default

    def pprint(self, val, length):
        lgd = ("%%-%ds" % length) % ("%s:" % self.desc[:length])
        sval = self.fmt % val
        return "%s%s" % (lgd, sval)


class Beam(ObjDebug):
    def copy(self, **nargs):
        t = self.__class__()
        t.__dict__.update(self.__dict__)
        t.__dict__.update(**nargs)
        return t

    def __init__(self, **nargs):
        for nnn, obj in self.__class__.__dict__.items():
            if isinstance(obj, Quantity):
                setattr(self, nnn, obj.default)

    def pprint(self, names, length=50):
        cls = self.__class__
        for nnn in names.split():
            val = getattr(self, nnn)
            try:
                lbl = getattr(cls, nnn).pprint(val, length=length)
            except AttributeError:
                lbl = ("%%-%ds%%s" % length) % ("%s:" % nnn, val)
            print(lbl)

    t_turnaround = 3
    cross_total = 110  # [mb]
    n_ip = 2
    lumi_lev = Quantity(5, "Leveling Luminosity [1e34 cm-2 s-1]")
    lumi_virt = Quantity(20, "Virtual  Luminosity [1e34 cm-2 s-1]")
    t_phys = Quantity(80, "Physics time [days]")
    n_bunches = Quantity(2736, "Number of Bunches")
    n_ppb = Quantity(2.2e11, "Protons per Bunch")
    lumi_ave_opt = Quantity(0, "Optimal Average Luminosity [1e34 cm-2 s-1]")
    lumi_int_opt = Quantity(0, "Optimal Integrated Luminosity [fb-1/year]")
    t_leveling = Quantity(10, "Leveling time [h]")
    t_decay = Quantity(3, "Luminosity decay time [h]")
    t_fill_max = Quantity(22, "Maximum fill time [h]")

    def mk_lumi_int_opt(self):
        """
        F. Zimmermann
        https://indico.cern.ch/getFile.py/access?contribId=4&sessionId=0&resId=1&materialId=slides&confId=175259
        """
        b = self.copy()
        teff = b.get_t_lifetime()
        cross_total = b.cross_total / 1e30
        ta = b.t_turnaround
        tmax = b.t_fill_max
        k = where(b.lumi_virt > b.lumi_lev, b.lumi_virt / b.lumi_lev, 1)
        llev = b.lumi_virt / k
        teff = b.n_bunches * b.n_ppb / b.n_ip / cross_total / llev / 3600 / 1e37
        osk = 1.0 / sqrt(k)
        tlvl = teff * (1 - osk)
        tmp_decay = -1 + osk + sqrt((1 - osk) ** 2 + (2 - osk) * ta / teff)
        tdec = teff / (2 - osk) * tmp_decay
        tlvl = where(tmax > tlvl, tlvl, tmax)
        tdec = where(tmax > tlvl + tdec, tdec, tmax - tlvl)
        tdec = where(tdec < 0, 0, tdec)
        tmp1 = tlvl + tdec * teff / (tdec + teff)
        tmp2 = tdec + tlvl + ta
        b.t_decay = tdec
        b.t_lifetime = teff
        b.t_leveling = tlvl
        b.t_fill_opt = tlvl + tdec
        b.lumi_ave_opt = llev * tmp1 / tmp2
        b.lumi_int_opt = b.lumi_ave_opt * b.t_phys * 24 * 3600 / 1e5
        return b

    def get_t_lifetime(b):
        cross_total = b.cross_total / 1e30
        teff = b.n_bunches * b.n_ppb / b.n_ip / cross_total / b.lumi_lev / 3600 / 1e37
        return teff

    def get_lumi_virt_from_t_leveling(b, tlev):
        return b.lumi_lev / (1 - tlev / b.get_t_lifetime()) ** 2

    def get_k(self):
        cross_total = b.cross_total / 1e30

    def mk_lumi_int_lev(self):
        b = self.copy()
        lfull = b.lumi_lev * b.t_phys * 24 * 3600 / 1e5  # [fb-1]
        b.lumi_int_lev = lfull * (
            self.t_leveling / (self.t_turnaround + self.t_leveling)
        )
        return b

    def plot_lumi_int_lev(self, cb=True):
        b = self.copy()
        lumi_lev = arange(0, 11.5, 0.5)
        t_leveling = arange(0, 11.1, 0.1)
        xx, yy = meshgrid(t_leveling, lumi_lev)
        b.t_leveling, b.lumi_lev = meshgrid(t_leveling, lumi_lev)
        b = b.mk_lumi_int_lev()
        zz = b.lumi_int_lev
        v = arange(0, 450, 50)
        p = plt.contourf(xx, yy, b.t_leveling, v, alpha=0.5)
        # plt.clabel(p,p.levels,inline=True,fmt="%d",fontsize=12)
        if cb:
            cbar = plt.colorbar()
            cbar.set_label(r"Integrated luminosity [$\rm fb^-1/year$]")
        plt.xlabel("Levelling duration [h]")
        plt.ylabel(r"Levelling Luminosity [$10^{34} \rm cm^{-2} s^{-1}$]")
        plt.title("no peak performance limits")
        plt.plot([0, 16], [5, 5], "k", lw=2)
        plt.plot([6, 6], [0, 11], "k", lw=2)
        xa, xb = plt.xlim(0, 11)
        ya, yb = plt.ylim(0, 11)
        plt.plot([xa, xb], [5, 5], "k", lw=2)
        plt.plot([6, 6], [ya, yb], "k", lw=2)
        return b

    def plot_lumi_int_opt(self, cb=True):
        b = self.copy()
        lumi_lev = arange(0, 10)
        lumi_virt = arange(0, 30, 0.1)
        xx, yy = meshgrid(lumi_virt, lumi_lev)
        b.lumi_virt, b.lumi_lev = xx, yy
        b = b.mk_lumi_int_opt()
        zz = b.lumi_int_opt
        v = arange(0, 450, 50)
        p = plt.contourf(b.t_fill_opt, yy, zz, v, alpha=0.5)
        # plt.clabel(p,p.levels,inline=True,fmt="%d",fontsize=12)
        if cb:
            cbar = plt.colorbar()
            cbar.set_label(r"Integrated luminosity [$\rm fb^-1/year$]")
        plt.xlabel("Fill duration [h]")
        plt.ylabel(r"Levelling Luminosity [$10^{34} \rm cm^{-2} s^{-1}$]")
        tlt = "Optimal luminosity: "
        tlt += r"$t_{\rm physics}$=%d d, " % b.t_phys
        tlt += r"$t_{\rm turnaround}$=%d h, " % b.t_turnaround
        tlt += r"$N_{\rm ppb}= %3.1f \cdot 10^{11}$" % (b.n_ppb / 1e11)
        plt.title(tlt)
        plt.xlim(1, 12)
        plt.ylim(1, 14)
        return b

    def plot_mlumi_int_opt(self, virt=[3, 5, 10, 15, 20, 22]):
        for vv in virt:
            self.plot_lumi_int_opt(vv, cb=False)

    def plot_lumi_int_virt(self, virt=20, n_ppb=2.2e11, n_bunches=2736, hl=5, cb=True):
        maxvirt = 11
        maxfill = 11
        b = self.copy()
        lumi_lev = linspace(0.001, min(maxvirt, virt), 30)
        b.lumi_virt = virt
        b.n_ppb = n_ppb
        b.n_bunches = n_bunches
        t_fill_max = linspace(0, maxfill, 30)
        xx, yy = meshgrid(t_fill_max, lumi_lev)
        b.t_fill_max, b.lumi_lev = xx, yy
        b = b.mk_lumi_int_opt()
        zz = b.lumi_int_opt
        v = arange(0, 450, 50)
        p = plt.contourf(b.t_fill_opt, yy, zz, v, alpha=0.5)
        # plt.clabel(p,p.levels,inline=True,fmt="%d",fontsize=12)
        if cb:
            cbar = plt.colorbar()
            cbar.set_label(r"Integrated luminosity [$\rm fb^{-1}/year$]")
        plt.xlabel("Fill duration [h]")
        plt.ylabel(r"Levelling Luminosity [$10^{34} \rm cm^{-2} s^{-1}$]")
        tlt = ""
        tlt += r"$L_{\rm virt}=%d  \cdot 10^{34}$, " % b.lumi_virt
        tlt += r"$N_{\rm ppb}= %3.1f \cdot 10^{11}$" % (b.n_ppb / 1e11)
        plt.title(tlt)
        xa, xb = plt.xlim(0, maxfill)
        ya, yb = plt.ylim(0, maxvirt)
        plt.plot([xa, xb], [hl, hl], "k", lw=2)
        # plt.plot([6,6],[ya,yb],'k',lw=2)
        plt.tight_layout()
        return b

    def plot_mlumi_int_opt(self, virt=[3, 5, 10, 15, 20, 22]):
        for vv in virt:
            self.plot_lumi_int_opt(vv, cb=False)

    def plot_lumi_int_emit(self, lumi_lev=5, turnaround=3, base=20, cb=True):
        b = self.copy()
        b.lumi_lev = lumi_lev
        b.t_turnaround = turnaround
        ppb = linspace(0.0, 3, 100) + 1e-6
        emit = linspace(0.7, 3, 100) + 1e-6
        xx, yy = meshgrid(ppb, emit)
        # emitf=1.1*yy+0.2*xx/yy + 0.10*xx**2/yy**0.5
        emitf = (
            1.1 * yy + 0.2 * xx / yy + 0.10 * xx**2 / yy**0.5
        )  # last term is fake
        b.lumi_virt = base * (xx * 0.95 / 2.2) ** 2 * (2.5 / emitf)
        b.n_ppb = xx * 1e11
        b = b.mk_lumi_int_opt()
        zz = b.lumi_int_opt
        # p=plt.contourf(xx,yy,b.t_lifetime,alpha=0.5)
        v = arange(0, 325, 25)
        p = plt.contourf(xx, yy, zz, v, alpha=0.5)
        # plt.clabel(p,p.levels,inline=True,fmt="%d",fontsize=12)
        if cb:
            cbar = plt.colorbar()
            cbar.set_label(r"Integrated luminosity [$\rm fb^{-1}/year$]")
        plt.xlabel("Injected bunch polulation [$10^{11}$ppb]")
        plt.ylabel(r"Injected emittance [$\mu$rad]")
        tlt = ""
        tlt += r"$L_{\rm level}=%3.2f  \cdot 10^{34}$, " % b.lumi_lev
        tlt += r"$k_{\rm IP1/5}= %4d$, " % b.n_bunches
        tlt += r"$t_{\rm ta}= %3.0f h$" % b.t_turnaround
        plt.title(tlt)
        # xa,xb=plt.xlim(0,3)
        # ya,yb=plt.ylim(0.5,3)
        # plt.plot([xa,xb],[5,5],'k',lw=2)
        # plt.plot([6,6],[ya,yb],'k',lw=2)
        return b

    def update_emit(self, x):
        pass

    # def __repr__(self):
    #  self.pprint("t_lifetime t_leveling t_decay")
    #  self.pprint("lumi_ave_opt lumi_int_opt")
