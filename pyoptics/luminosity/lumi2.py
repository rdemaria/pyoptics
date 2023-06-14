import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, e, m_p, epsilon_0
from scipy.special import erf

import scipy.integrate
import scipy.optimize

clight = 299792458
pmass = 0.938272013e9


def ftoint(u, alpha=1.0, Ax=1.0, Ay=1.0):
    """alpha=theta_c/2 * sqrt(emit/betxstar)
    Ax=sigma_z/betxstar
    Ay=sigma_z/betystar
    """
    dx = 1 + (Ax * u) ** 2
    dy = 1 + (Ay * u) ** 2
    n = np.exp(-(u**2) * (1 + (alpha * Ax) ** 2 / dx))
    d = np.sqrt(dx * dy)
    return n / d


def mkint(alpha=290e-6, Ax=1.0, Ay=1.0, debug=False):
    i = scipy.integrate.quad(ftoint, 0, np.inf, args=(alpha, Ax, Ay))
    if debug:
        print("Integral tolerance: %e" % i[1])
    return 2 / np.sqrt(np.pi) * i[0]


def gaussian_long_profile(s, ct, sigz, thetax, ct01, ct02):
    costh = np.cos(thetax * 0.5)
    return np.exp(
        -((s * costh - (ct + ct01)) ** 2 + (s * costh + (ct + ct02)) ** 2)
        / (2 * sigz**2)
    )


def Fxy(
    s,
    ct,
    betx,
    bety,
    emitx,
    emity,
    theta,
    ccx,
    ccy,
    ctx_p,
    ctx_m,
    cty_p,
    cty_m,
    kx,
    ky,
):
    """
    s : longitudinal coordinate
    ct : integration variable
    betx : beta function in the crossing plane
    bety : beta function in the separation plane
    emitx : emittance in the crossing plane
    emity : emittance in the separation plane
    ccx: crab cavity angle in the crossing plane
    ccy: crab cavity angle in the separation plane
    theta: crossing angle
    """
    sigx = np.sqrt(betx * emitx)
    sigy = np.sqrt(bety * emity)
    bsx = 1 + (s / betx) ** 2  # betx(s)=bsx * betx
    bsy = 1 + (s / bety) ** 2  # bety(s)=bsy * bety
    bsx_sq = np.sqrt(bsx)
    bsy_sq = np.sqrt(bsy)

    argx = np.sin(kx * (s + ctx_m / 2)) * np.cos(kx * (ct + ctx_p / 2)) * np.sin(
        ccx / 2
    ) / kx - s * np.sin(theta / 2)
    argy = (
        np.cos(ky * (s + cty_m / 2))
        * np.sin(ky * (ct + cty_p / 2))
        * np.sin(ccy / 2)
        / ky
    )

    res = (
        np.exp(-((argx / sigx) ** 2) / bsx)
        / bsx_sq
        * np.exp(-((argy / sigy) ** 2) / bsy)
        / bsy_sq
    )
    return res


class GaussianLongDist:
    def __init__(self, sigz):
        self.sigz = sigz


class QGaussianLongDist:
    def __init__(self, sigz, q):
        self.sigz = sigz
        self.q = q


class InteractionPoint:
    """Interaction Point class
    name: name of the IP
    betx [m]: beta function at the IP in the horizontal plane
    bety [m]: beta function at the IP in the vertical plane
    sepx [m]: horizontal separation at the IP
    sepy [m]: vertical separation at the IP
    thetax [rad]: horizontal crossing angle at the IP
    thetay [rad]: vertical crossing angle at the IP
    dx[m]: horizontal dispersion at the IP
    dy[m]: vertical dispersion at the IP
    ccx[rad]: horizontal crab cavity angle at the IP
    ccy[rad]: vertical crab cavity angle at the IP
    ccrf[rad]: crab cavity RF frequency at the IP
    cclag[rad]: crab cavity RF phase lag at the IP
    visible_cross_section[mb]: visible cross section at the IP
    total_cross_section[mb]: total cross section at the IP
    """

    def __init__(
        self,
        name="ip5",
        betx=1,
        bety=1,
        sepx=0,
        sepy=0,
        thetax=0,
        thetay=0,
        dx=0,
        dy=0,
        ccx=0,
        ccy=0,
        ccrf=0,
        cclag=0,
        visible_cross_section=81,
        total_cross_section=81,
    ):
        self.name = name
        self.betx = betx
        self.bety = bety
        self.sepx = sepx
        self.sepy = sepy
        self.thetax = thetax
        self.thetay = thetay
        self.dx = dx
        self.dy = dy
        self.ccx = ccx
        self.ccy = ccy
        self.ccrf = ccrf
        self.cclag = cclag
        self.visible_cross_section = visible_cross_section
        self.total_cross_section = total_cross_section

    def normalized_crossing_angle(self, bunch):
        """Normalized crossing angle"""
        phix = (self.thetax - self.ccx) / (np.sqrt(bunch.emitx / self.betx))
        phiy = (self.thetay - self.ccy) / (np.sqrt(bunch.emity / self.bety))
        return phix, phiy

    def geometric_factor(self, bunch):
        """Geometric factor"""
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        effx = self.thetax - self.ccx
        effy = self.thetay - self.ccy
        geox = 1 / np.sqrt(1 + ((bunch.sigz * effx) / (2 * sigx)) ** 2)
        geoy = 1 / np.sqrt(1 + ((bunch.sigz * effy) / (2 * sigy)) ** 2)
        return geox, geoy

    def lumi_headon(self, bunch):
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        L0 = (bunch.ppb**2 * bunch.nb * bunch.frev) / (4 * np.pi * sigx * sigy)
        return L0

    def lumi(self, bunch):
        if self.ccrf == 0:
            if min(self.betx,self.bety) > 2*bunch.sigz:
                return self.lumi_simple(bunch)
            else:
                return self.lumi_hourglass_noccrf(bunch)
        else:
            return self.lumi_ccrf(bunch)
        
    def lumi_simple(self, bunch):
        """Luminosity at the IP"""
        L0 = self.lumi_headon(bunch)
        effx = self.thetax - self.ccx
        effy = self.thetay - self.ccy
        theta = np.sqrt(effx**2 + effy**2)
        xp = np.arctan2(effx, effy)
        sigx= np.sqrt(self.betx * bunch.emitx)
        #assert xp >= 0
        beta_cross = abs(self.betx * np.cos(xp) + self.bety * np.sin(xp))
        beta_sep = abs(self.betx * np.sin(xp) + self.bety * np.cos(xp))
        emit_cross = abs(bunch.emitx * np.cos(xp) + bunch.emity * np.sin(xp))
        emit_sep = abs(bunch.emitx * np.sin(xp) + bunch.emity * np.cos(xp))
        sig_cross = np.sqrt(beta_cross * emit_cross)
        sep_cross = self.sepx * np.cos(xp) + self.sepy * np.sin(xp)
        sep_sep = self.sepx * np.sin(xp) + self.sepy * np.cos(xp)
        print(sep_sep,sep_cross)
        sig_sep = np.sqrt(beta_sep * emit_sep)
        sig_cross = np.sqrt(beta_cross * emit_cross)
        print(sig_sep,sig_cross)
        geo = 1 / np.sqrt(1 + ((bunch.sigz * theta) / (2 * sig_cross)) ** 2)
        #w=np.exp(-(sep_cross/2/sig_cross)**2)*np.exp(-(sep_sep/2/sig_sep)**2)
        print(self.sepx,sigx,self.sepx/2/sigx)
        w=np.exp(-(self.sepx/2/sigx)**2)
        a=(np.sin(theta/2)/sig_cross)**2 + (np.cos(theta/2)/bunch.sigz)**2
        b=sep_cross*np.sin(theta/2)/2/sig_cross**2
        print(L0,geo,w,a,b)
        return L0* geo * w * np.exp(b**2/a)

        alpha = theta / (2 * np.sqrt(emit / beta_cross))
        Ax = bunch.sigz / beta_cross
        Ay = bunch.sigz / beta_sep
        factor = mkint(alpha, Ax, Ay)
        LL = L0 * factor
        return LL
    

    def lumi_hourglass_noccrf(self, bunch):
        """Luminosity at the IP"""
        L0 = self.lumi_headon(bunch)
        effx = self.thetax - self.ccx
        effy = self.thetay - self.ccy
        theta = np.sqrt(effx**2 + effy**2)
        xp = np.arctan2(effx, effy)
        assert xp >= 0
        beta_cross = self.betx * np.cos(xp) + self.bety * np.sin(xp)
        beta_sep = self.betx * np.sin(xp) + self.bety * np.cos(xp)
        emit = bunch.emitx * np.cos(xp) + bunch.emity * np.sin(xp)
        alpha = theta / (2 * np.sqrt(emit / beta_cross))
        Ax = bunch.sigz / beta_cross
        Ay = bunch.sigz / beta_sep
        factor = mkint(alpha, Ax, Ay)
        LL = L0 * factor
        return LL

    def lumi_ccrf(self, bunch, debug=False, long_profile=gaussian_long_profile):
        """
        NB: don't work with tilted crossing and coupling
        """

        # pseudo rotation of the crossing plane
        theta = np.sqrt(self.thetax**2 + self.thetay**2)  # crossing angle
        phi = np.arctan2(abs(self.thetax), abs(self.thetay))  # crossing plane
        cc = np.cos(phi)
        ss = np.sin(phi)
        assert np.close(phi % (np.pi / 2), 0)
        beta_cross = self.betx * cc + self.bety * ss
        beta_sep = self.betx * ss + self.bety * cc
        cc_cross = self.ccx * cc + self.ccy * ss
        cc_sep = self.ccx * ss + self.ccy * cc
        emit_cross = bunch.emitx * cc + bunch.emity * ss
        emit_sep = bunch.emitx * ss + bunch.emity * cc

        kx = self.ccrf / clight * 2 * np.pi
        ky = self.ccrf / clight * 2 * np.pi

        ctx_p = 0  # self.cctx[1] + self.cctx[0]
        ctx_m = 0  # self.cctx[1] - self.cctx[0]

        cty_p = 0  # self.ccty[1] + self.ccty[0]
        cty_m = 0  # self.ccty[1] - self.ccty[0]

        ftoint = lambda s, ct: Fxy(
            s,
            ct,
            beta_cross,
            beta_sep,
            emit_cross,
            emit_sep,
            theta,
            cc_cross,
            cc_sep,
            ctx_p,
            ctx_m,
            cty_p,
            cty_m,
            kx,
            ky,
        ) * long_profile(s, bunch.sigz)

        res = scipy.integrate.dblquad(
            ftoint,
            -np.inf,
            np.inf,
            -np.inf,
            np.inf,
        )

        if debug:
            print("Integral tolerance: %e" % res[1])

        LO = self.lumi_headon(bunch)
        return res[0] * LO * (np.cos(theta / 2)) / (np.pi * bunch.sigz**2)


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
        ppb=2.3e11,
        emitnx=2.5e-6,
        emitny=2.5e-6,
        sigz=7.61e-2,
        sigdpp=1.2e-4,
        energy=7e12,
        ips=(),
        longdist=GaussianLongDist,
        frev=11245.5,
        pmass=0.9382720813e9,
    ):
        self.nb = nb
        self.ppb = ppb
        self.emitnx = emitnx
        self.emitny = emitny
        self.sigz = sigz
        self.sigdpp = sigdpp
        self.energy = energy
        self.ips = ips
        self.longdist = longdist
        self.frev = frev
        self.pmass = pmass

    @property
    def gamma(self):
        return self.energy / self.pmass

    @property
    def emitx(self):
        return self.emitnx / self.gamma

    @property
    def emity(self):
        return self.emitny / self.gamma

    def lumi(self):
        return np.array([ip.lumi(self) for ip in self.ips])


class LevellingProcess:
    def __init__(self, bunches, interactions, lumi_start, lumi_lev, lumi_ramp_time):
        self.bunches = bunches
        self.interactions = interactions
        self.lumi_start = lumi_start
        self.lumi_lev = lumi_lev
        self.lumi_ramp_time = lumi_ramp_time


class LumiCalculator:
    proton_mass = 0.938  # in GeV
    emitt = 2.5e-6  # [m]
    frev = 299792458 / 26658.8832  # hz
    t1 = 0  # -50e-12# bunch time delay [s]
    t2 = 0  # 10e-12# bunch time delay [s]
    # wcc = 400.0 * 1e6 # RF frequency ???

    def __init__(
        self,
        crossing_angle,
        crab_angle,
        energy,
        nbunches,
        nparticles,
        betx,
        bety,
        sig_s,
        wcc=400 * 1e6,
    ):
        self.wcc = wcc
        self.crab_angle = crab_angle  # in radians
        self.crossing_angle = crossing_angle  # rad
        self.Nb = nbunches
        self.Np = nparticles
        self.gamma = energy / self.proton_mass
        self.geom_emitt = self.emitt / self.gamma

        self.betx = betx
        self.bety = bety

        self.sigx = np.sqrt(self.geom_emitt * self.betx)
        self.sigy = np.sqrt(self.geom_emitt * self.bety)
        self.sig_s = sig_s
