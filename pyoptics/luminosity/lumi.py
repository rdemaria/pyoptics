import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c as clight

import scipy.optimize

from . import lumi_guido

twopi = 2 * np.pi
fourpi = 4 * np.pi


def crossing_spec(x1=0, x2=0, y1=0, y2=0, px1=0, px2=0, py1=0, py2=0, degree=False):
    xingh = px1 - px2
    xingv = py1 - py2
    seph = x1 - x2
    sepv = y1 - y2
    angh = px1 + px2
    angv = py1 + py2
    offh = x1 + x2
    offv = y1 + y2
    phi = np.arctan2(xingv, xingh)
    xing = np.sqrt(xingh**2 + xingv**2)
    sep = np.sqrt(seph**2 + sepv**2)
    ang = np.sqrt(angh**2 + angv**2)
    off = np.sqrt(offh**2 + offv**2)
    if degree:
        phi *= 180 / np.pi
    return xing / 2, sep / 2, ang / 2, off / 2, phi


def orbit_spec(xing=0, sep=0, ang=0, off=0, phi=0, degree=False):
    if degree:
        phi *= np.pi / 180
    cc = np.cos(phi)
    ss = np.sin(phi)
    xingh = cc * xing
    xingv = ss * xing
    offh = -ss * off
    offv = cc * off
    seph = -ss * sep
    sepv = cc * sep
    angh = -ss * ang
    angv = cc * ang
    x1 = seph + offh
    x2 = -seph + offh
    y1 = sepv + offv
    y2 = -sepv + offv
    px1 = xingh + angh
    px2 = -xingh + angh
    py1 = xingv + angv
    py2 = -xingv + angv
    return x1, x2, y1, y2, px1, px2, py1, py2


class IP:
    """Interaction Point class
    name: name of the IP
    betx [m]: beta function at the IP in the horizontal plane
    bety [m]: beta function at the IP in the vertical plane
    sepx [m]: half horizontal separation at the IP
    sepy [m]: half vertical separation at the IP
    px [rad]: half horizontal crossing angle at the IP
    py [rad]: half vertical crossing angle at the IP
    dx[m]: half horizontal dispersion at the IP
    dy[m]: half vertical dispersion at the IP
    ccvx[rad]: half horizontal crab cavity voltage at the IP
    ccvy[rad]: half vertical crab cavity voltage at the IP
    ccrf[rad]: crab cavity RF frequency at the IP
    cclag[rad]: crab cavity RF phase lag at the IP
    visible_cross_section[mb]: visible cross section at the IP
    total_cross_section[mb]: total cross section at the IP
    """

    def set_crossing(self, xing, sep, ang, off, phi, degree=True):
        (
            self.x1,
            self.x2,
            self.y1,
            self.y2,
            self.px1,
            self.px2,
            self.py1,
            self.py2,
        ) = orbit_spec(xing, sep, ang, off, phi, degree=degree)

    def get_crossing(self, degree=True):
        return crossing_spec(
            self.x1,
            self.x2,
            self.y1,
            self.y2,
            self.px1,
            self.px2,
            self.py1,
            self.py2,
            degree=degree,
        )

    def get_phi(self, degree=True):
        phi = np.arctan2(self.py, self.px)
        if degree:
            phi *= 180 / np.pi
        return phi

    @property
    def phi(self):
        return self.get_phi()

    @phi.setter
    def phi(self, phi):
        self.set_crossing(*self.get_crossing(degree=True), phi=phi, degree=True)

    @property
    def sep(self):
        return np.sqrt(self.sepx**2 + self.sepy**2)

    def sep_(self, sep):
        phi_rad = self.phi / 180 * np.pi
        cc = np.cos(phi_rad)
        ss = np.sin(phi_rad)
        self.sepx = -ss * sep
        self.sepy = cc * sep
        return self

    sep = sep.setter(sep_)

    @property
    def sigx(self):
        return np.sqrt(self.betx * self.emitx + self.dx**2)

    @property
    def sigy(self):
        return np.sqrt(self.bety * self.emity + self.dy**2)

    def __init__(
        self,
        name="ip5",
        betx=1,
        bety=1,
        alfx=0,
        alfy=0,
        sepx=0,
        sepy=0,
        px=0,
        py=0,
        dx=0,
        dy=0,
        ccx=0,
        ccy=0,
        # ccr12=23.3305,
        # ccr34=23.3305,
        ccrf=400e6,
        cclag=0,
        visible_cross_section=81,
        total_cross_section=100,
    ):
        self.name = name
        self.betx = betx
        self.bety = bety
        self.alfx = alfx
        self.alfy = alfy
        self.sepx = sepx
        self.sepy = sepy
        self.px = px
        self.py = py
        self.dx = dx
        self.dy = dy
        self.dpx = 0
        self.dpy = 0
        self.ccx = ccx
        self.ccy = ccy
        # self.ccr12 = ccr12
        # self.ccr34 = ccr34
        self.ccrf = ccrf
        self.cclag = cclag
        self.visible_cross_section = visible_cross_section
        self.total_cross_section = total_cross_section

    def __repr__(self) -> str:
        out = []
        for kk in ["name", "betx", "bety", "sepx", "sepy", "px", "py", "ccx", "ccy"]:
            vv = getattr(self, kk)
            if vv != 0:
                out.append(f"{kk}={vv!r}")
        return f"IP({', '.join(out)})"

    def clone(self, **kwargs):
        ip = IP(
            name=self.name,
            betx=self.betx,
            bety=self.bety,
            alfx=self.alfx,
            alfy=self.alfy,
            sepx=self.sepx,
            sepy=self.sepy,
            px=self.px,
            py=self.py,
            dx=self.dx,
            dy=self.dy,
            ccx=self.ccx,
            ccy=self.ccy,
            # ccr12=self.ccr12,
            # ccr34=self.ccr34,
            ccrf=self.ccrf,
            cclag=self.cclag,
            visible_cross_section=self.visible_cross_section,
            total_cross_section=self.total_cross_section,
        )
        for k, v in kwargs.items():
            setattr(ip, k, v)
        return ip

    def pileup(self, bunch):
        """Pile-up"""
        l = self.luminosity(bunch)
        return l * self.visible_cross_section * 1e-31 / (bunch.nb * bunch.frev)

    def normalized_crossing_angle(self, bunch):
        """Normalized crossing angle"""
        phix = (self.px) / (np.sqrt(bunch.emitx / self.betx))
        phiy = (self.py) / (np.sqrt(bunch.emity / self.bety))
        return phix, phiy

    def normalized_separation(self, bunch):
        """Normalized separation"""
        nsepx = self.sepx / (np.sqrt(bunch.emitx * self.betx))
        nsepy = self.sepy / (np.sqrt(bunch.emity * self.bety))
        return nsepx, nsepy

    # def crabing_angles(self, bunch):
    #     """Crabbing angles"""
    #     phix = self.ccr12 * self.ccvx / bunch.energy * twopi * self.ccrf / clight
    #     phiy = self.ccr34 * self.ccvy / bunch.energy * twopi * self.ccrf / clight
    #     return phix, phiy

    def geometric_factor(self, bunch):
        """Geometric factor"""
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        effx = self.px + self.ccx
        effy = self.py + self.ccy
        geox = 1 / np.sqrt(1 + ((bunch.sigz * effx) / sigx) ** 2)
        geoy = 1 / np.sqrt(1 + ((bunch.sigz * effy) / sigy) ** 2)
        return geox, geoy

    def separation_factor(self, bunch):
        """Separation factor"""
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        fx = np.exp(-self.sepx**2 / sigx**2)
        fy = np.exp(-self.sepy**2 / sigy**2)
        return fx, fy

    def lumi_headon(self, bunch):
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        L0 = (bunch.ppb**2 * bunch.nb * bunch.frev) / (fourpi * sigx * sigy)
        return L0

    def lumi_simple(self, bunch):
        L0 = self.lumi_headon(bunch)
        geox, geoy = self.geometric_factor(bunch)
        fx, fy = self.separation_factor(bunch)
        L = L0 * geox * geoy * fx * fy
        return L

    def luminosity(self, bunch, verbose=False):
        ccr12 = 1
        ccr34 = 1
        ccvx = self.ccx / ccr12 * bunch.energy / twopi / self.ccrf * clight
        ccvy = self.ccy / ccr34 * bunch.energy / twopi / self.ccrf * clight

        return lumi_guido.luminosity(
            f=bunch.frev,
            nb=bunch.nb,
            N1=bunch.ppb,
            N2=bunch.ppb,
            x_1=self.sepx,
            x_2=-self.sepx,
            y_1=self.sepy,
            y_2=-self.sepy,
            px_1=self.px,
            px_2=-self.px,
            py_1=self.py,
            py_2=-self.py,
            energy_tot1=bunch.energy,
            energy_tot2=bunch.energy,
            deltap_p0_1=bunch.delta,
            deltap_p0_2=bunch.delta,
            epsilon_x1=bunch.emitnx,
            epsilon_x2=bunch.emitnx,
            epsilon_y1=bunch.emitny,
            epsilon_y2=bunch.emitny,
            sigma_z1=bunch.sigz,
            sigma_z2=bunch.sigz,
            beta_x1=self.betx,
            beta_x2=self.betx,
            beta_y1=self.bety,
            beta_y2=self.bety,
            alpha_x1=self.alfx,
            alpha_x2=-self.alfx,
            alpha_y1=self.alfy,
            alpha_y2=-self.alfy,
            dx_1=self.dx,
            dx_2=self.dx,
            dy_1=self.dy,
            dy_2=self.dy,
            dpx_1=self.dpx,
            dpx_2=self.dpx,
            dpy_1=self.dpy,
            dpy_2=self.dpy,
            CC_V_x_1=ccvx,
            CC_f_x_1=self.ccrf,
            CC_phase_x_1=0,
            CC_V_x_2=-ccvx,
            CC_f_x_2=self.ccrf,
            CC_phase_x_2=0,
            CC_V_y_1=ccvy,
            CC_f_y_1=self.ccrf,
            CC_phase_y_1=0,
            CC_V_y_2=-ccvy,
            CC_f_y_2=self.ccrf,
            CC_phase_y_2=0,
            R12_1=ccr12,
            R22_1=0,
            R34_1=ccr34,
            R44_1=0,
            R12_2=ccr12,
            R22_2=0,
            R34_2=ccr34,
            R44_2=0,
            verbose=verbose,
            sigma_integration=3,
        )

    def betastar_from_lumi(self, bunch, target, betaratio=1):
        """Solve for the betastar that give a target luminosity"""
        ip = self.clone()

        def ftosolve(beta):
            ip.betx = beta
            ip.bety = beta * betaratio
            return ip.luminosity(bunch) - target

        res = scipy.optimize.root(ftosolve, ip.betx)
        if res.success:
            ftosolve(res.x[0])
            return ip
        else:
            print(res)
            raise ValueError("Could not find betastar")

    def betastar_from_pileup(self, bunch, target, betaratio=1):
        """Solve for the betastar that give a target pileup"""
        ip = self.clone()

        def ftosolve(beta):
            ip.betx = beta
            ip.bety = beta * betaratio
            return ip.pileup(bunch) - target

        res = scipy.optimize.root(ftosolve, ip.betx)
        if res.success:
            beta = res.x[0]
            ip.betx = beta
            ip.bety = beta * betaratio
            return ip
        else:
            print(res)
            raise ValueError("Could not find betastar")

    def sep_from_lumi(self, bunch, target):
        ip = self.clone()

        def ftosolve(sep):
            ip.sep_(sep)
            return ip.luminosity(bunch) - target

        res = scipy.optimize.root(ftosolve, ip.sep)
        if res.success:
            ip.sep_(res.x[0])
            return ip
        else:
            print(res)
            raise ValueError("Could not find separation")

    def sep_from_lumi_simple(self, bunch, target):
        l_nosep = self.clone().sep_(0).lumi_simple(bunch)
        factor = target / l_nosep
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        if factor > 1:
            sep = 0  # no levelling possible
        else:
            sep = np.sqrt(-np.log(factor) * sigx * sigy)  # small error
        return self.clone().sep_(sep)

    def burnoff(self, bunch):
        """Proton per bunch loss in 1 sec"""
        return self.luminosity(bunch) * self.total_cross_section * 1e-31 / bunch.nb

    def info(self, bunch):
        print(f"sigma_x                  : {np.sqrt(self.betx * bunch.emitx)}")
        print(f"sigma_y                  : {np.sqrt(self.bety * bunch.emity)}")
        print(f"Normalized crosing angles: {self.normalized_crossing_angle(bunch)}")
        print(f"Normalized separations   : {self.normalized_separation(bunch)}")


class Bunch:
    """Bunch class

    nb: number of particles per bunch
    ppb: number of protons per bunch
    emitx [m.rad]: horizontal emittance
    emity [m.rad]: vertical emittance
    sigz [m]: RMS bunch length
    sigdpp: RMS energy spread
    energy [eV]: beam energy
    ips: list of InteractionPoint objects
    longdist: longitudinal distribution
    pmass [eV]: proton mass
    frev [Hz]: revolution frequency
    """

    def __init__(
        self,
        nb=1960,
        ppb=2.2e11,
        emitnx=2.5e-6,
        emitny=2.5e-6,
        sigz=7.61e-2,
        sigdpp=1.2e-4,
        energy=7e12,
        ips=(),
        long_dist="gaussian",
        frev=11245.5,
        pmass=0.9382720813e9,
        delta=1,
    ):
        self.nb = nb
        self.ppb = ppb
        self.emitnx = emitnx
        self.emitny = emitny
        self.sigz = sigz
        self.sigdpp = sigdpp
        self.energy = energy
        self.ips = ips
        self.longdist = long_dist
        self.frev = frev
        self.pmass = pmass
        self.delta = delta

    def __repr__(self) -> str:
        out = []
        for kk in ["energy", "nb", "ppb", "emitnx", "emitny", "sigz"]:
            vv = getattr(self, kk)
            if vv != 0:
                out.append(f"{kk}={vv!r}")
        return f"IP({', '.join(out)})"

    @property
    def gamma(self):
        return self.energy / self.pmass

    @property
    def emitx(self):
        return self.emitnx / self.gamma

    @property
    def emity(self):
        return self.emitny / self.gamma

    def luminosity(self, verbose=False):
        return np.array([ip.luminosity(self, verbose=verbose) for ip in self.ips])

    def ip_info(self):
        for ip in self.ips:
            ip.info(self)

    def clone(self, **kwargs):
        bb = Bunch(
            nb=self.nb,
            ppb=self.ppb,
            emitnx=self.emitnx,
            emitny=self.emitny,
            sigz=self.sigz,
            sigdpp=self.sigdpp,
            energy=self.energy,
            ips=self.ips,
            long_dist=self.longdist,
            frev=self.frev,
            pmass=self.pmass,
            delta=self.delta,
        )
        for k, v in kwargs.items():
            setattr(bb, k, v)
        return bb


class BetaStarLeveling:
    def __init__(
        self,
        bunches,
        ips,
        lumi_start,
        lumi_lev,
        lumi_ramp_time,
        betastar_end=0.15,
        lumi2=1.4e35,
        lumi8=2e37,
    ):
        self.bunches = bunches
        self.ips = ips
        self.lumi_start = lumi_start
        self.lumi_lev = lumi_lev
        self.lumi_ramp_time = lumi_ramp_time
        self.betastar_end = betastar_end
        self.lumi2 = lumi2
        self.lumi8 = lumi8

    def get_colliding_bunches(self, verbose=True):
        nc = {1: 0, 2: 0, 5: 0, 8: 0}
        for ip in (1, 2, 5, 8):
            for bb in self.bunches:
                if ip in bb.ips:
                    nc[ip] += bb.nb
        return nc

        ips = (ip1, ip2, ip5, ip8)
        if verbose:
            print(f"{ip1.name:6} beta {ip1.betx:5.2f},{ip1.bety:5.2f} m")
            print(f"{ip5.name:6} beta {ip5.betx:5.2f},{ip5.bety:5.2f} m")
            print(f"{ip2.name:6} sep  {ip2.sepx*1e3:5.2f} mm")
            print(f"{ip8.name:6} sep  {ip2.sepy*1e3:5.2f} mm")
        return ips

    def level_betasep(self, ips, lumis, bunch, verbose=True):
        ip1, ip5, ip2, ip8 = ips
        bunch = self.bunches[0]
        nc = self.get_colliding_bunches(verbose=verbose)
        lumi1 = lumis[0] * bunch.nb / nc[1]
        lumi5 = lumis[1] * bunch.nb / nc[5]
        lumi2 = lumis[2] * bunch.nb / nc[2]
        lumi8 = lumis[3] * bunch.nb / nc[8]

        ip1 = self.ips[1].betastar_from_lumi(bunch, lumi1)
        ip5 = self.ips[5].betastar_from_lumi(bunch, lumi5)
        ip2 = self.ips[2].sep_from_lumi_simple(bunch, lumi2)
        ip8 = self.ips[8].sep_from_lumi_simple(bunch, lumi8)

        if verbose:
            print(ip1.name, ip1.betx, ip1.bety, ip1.luminosity(bunch))
            print(ip5.name, ip5.betx, ip5.bety, ip5.luminosity(bunch))
            print(ip2.name, ip2.sepx, ip2.luminosity(bunch))
            print(ip8.name, ip8.sepy, ip8.luminosity(bunch))

        return ip1, ip5, ip2, ip8

    def betastar_leveling(self, fillduration=15 * 3600, dt=60, verbose=True):
        """Integrate the luminosity over time"""
        out = []

        # find initial conditions
        ip1, ip5, ip2, ip8 = self.level_betasep(
            self.ips,
            (self.lumi_start, self.lumi_start, self.lumi2, self.lumi8),
            self.bunches[0],
            verbose=True,
        )

        # first step
        bunches = [b.clone() for b in self.bunches]
        tt = 0  # luminosity time
        ips = {1: ip1, 2: ip2, 5: ip5, 8: ip8}
        out.append((tt, bunches, ips))
        # lumi ramp
        lumi15 = self.lumi_start
        dlumi = (self.lumi_lev - self.lumi_start) / self.lumi_ramp_time * dt
        dcc = 190e-6 / self.lumi_ramp_time * dt
        while lumi15 < self.lumi_lev:
            if verbose:
                print(out[-1])
            # clone and burnoff
            bunches = [b.clone() for b in bunches]
            ips = {1: ip1.clone(), 2: ip2.clone(), 5: ip5.clone(), 8: ip8.clone()}
            for bb in bunches:
                for ip in bb.ips:
                    bb.ppb -= ips[ip].burnoff(bb) * dt
                if verbose:
                    print(bb.ppb)
            # levelling
            lumi15 += dlumi
            ips[1].ccy += dcc
            ips[5].ccx += dcc
            ip1, ip5, ip2, ip8 = self.level_betasep(
                ips,
                (lumi15, lumi15, self.lumi2, self.lumi8),
                self.bunches[0],
                verbose=True,
            )

            # update
            out.append((tt, bunches, ip1, ip5, ip2, ip8))
            tt += dt
            if tt > fillduration:
                break
        # lumi leveling
        lumi15 = self.lumi_lev  # small error here
        while ip1.betx > self.betastar_end:
            if verbose:
                print(out[-1])
            # clone and burnoff
            bunches = [b.clone() for b in bunches]
            ips = {1: ip1.clone(), 2: ip2.clone(), 5: ip5.clone(), 8: ip8.clone()}
            for bb in bunches:
                for ip in bb.ips:
                    bb.ppb -= ips[ip].burnoff(bb) * dt
            # levelling
            ip1, ip5, ip2, ip8 = self.level_betasep(
                self.ips, (lumi15, lumi15, self.lumi2, self.lumi8), self.bunches[0]
            )
            # update
            out.append((tt, bunches, ip1, ip5, ip2, ip8))
            tt += dt
            if tt > fillduration:
                break
        # lumi decay
        while tt < fillduration:
            if verbose:
                print(out[-1])
            # clone and burnoff
            bunches = [b.clone() for b in bunches]
            ips = {1: ip1.clone(), 2: ip2.clone(), 5: ip5.clone(), 8: ip8.clone()}
            for bb in bunches:
                for ip in bb.ips:
                    bb.ppb -= ips[ip].burnoff(bb) * dt
            # levelling
            ip1, ip5, ip2, ip8 = self.level_betasep(
                self.ips,
                (lumi15, lumi15, self.lumi2, self.lumi8),
                self.bunches[0],
                verbose=True,
            )
            # update
            out.append((tt, bunches, ip1, ip5, ip2, ip8))
            tt += dt
        return StableBeam(*zip(*out))


class StableBeam:
    def __init__(self, tt, bunches, ip1, ip5, ip2, ip8):
        self.tt = tt
        self.bunches = bunches
        self.ip1 = ip1
        self.ip5 = ip5
        self.ip2 = ip2
        self.ip8 = ip8
