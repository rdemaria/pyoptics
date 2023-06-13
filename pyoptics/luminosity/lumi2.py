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
        """Luminosity at the IP"""
        sigx = np.sqrt(self.betx * bunch.emitx)
        sigy = np.sqrt(self.bety * bunch.emity)
        L0 = (bunch.ppb**2 * bunch.nb * bunch.frev) / (4 * np.pi * sigx * sigy)
        effx = self.thetax - self.ccx
        effy = self.thetay - self.ccy
        theta=np.sqrt(effx**2+effy**2)
        xp=np.arctan2(effx,effy)
        assert xp>=0
        beta_cross=self.betx*np.cos(xp)+self.bety*np.sin(xp)
        beta_sep=self.betx*np.sin(xp)+self.bety*np.cos(xp)
        emit=bunch.emitx*np.cos(xp)+bunch.emity*np.sin(xp)
        alpha = theta / (2 * np.sqrt(emit / beta_cross))
        Ax = bunch.sigz / beta_cross
        Ay = bunch.sigz / beta_sep
        factor = mkint(alpha, Ax, Ay)
        LL = L0 * factor
        return LL


    def L0(self):
        return self.frev*self.Nb*self.Np*self.Np/(4*np.pi*self.sigx*self.sigy)*1e-4


    def Fxy(self, s, ct, ctx_p, ctx_m, cty_p, cty_m, kx, rfc=True):
        sigmax2 = self.sigx**2
        sigmay2 = self.sigy**2
        bsx = (1+(s/self.betx)**2)
        bsy = (1+(s/self.bety)**2)
        bsx_sq = np.sqrt(bsx)
        bsy_sq = np.sqrt(bsy)

        if rfc is True:
            res = np.exp( -(np.sin(kx*(s+ctx_m*0.5))*np.cos(kx*(ct+ctx_p*0.5))*np.sin(self.crab_angle*0.5)/kx -s*np.sin(self.crossing_angle*0.5))**2 / (sigmax2 * bsx) ) / bsx_sq / bsy_sq
        elif rfc is False:        
            res = np.exp( -(s*np.sin(self.crab_angle*0.5) - s*np.sin(self.crossing_angle*0.5))**2 / (sigmax2 * bsx) ) / bsx_sq / bsy_sq   
        return res


    def long_profile(self, s, ct, ct01, ct02):
        costh = np.cos(self.crossing_angle*0.5)
        return np.exp( -((s*costh - (ct + ct01))**2 + (s*costh + (ct + ct02))**2) / (2*self.sig_s**2) )
 
    def luminosity_density(self, s, ct, ctx=(0,0), cty=(0,0), rfc=True):
        """
        s : longitudinal coordinate
        ct : integration variable
        ctx : time lag of the crab cavity in the x plane
        cty : time lag of the crab cavity in the y plane
        """

        kx = self.wcc/clight*2*np.pi

        ctx_p = ctx[1] + ctx[0]
        ctx_m = ctx[1] - ctx[0]

        cty_p = cty[1] + cty[0]
        cty_m = cty[1] - cty[0]

        fxy = self.Fxy(s, ct, ctx_p, ctx_m, cty_p, cty_m, kx, rfc)

        integrand = fxy*self.long_profile(s, ct, 30e-12*clight, 30e-12*clight)

        return integrand * (self.L0()*np.cos(self.crossing_angle*0.5))/ (np.pi*self.sig_s**2)

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


class LevellingProcess:
    def __init__(self, bunches, interactions, lumi_start, lumi_lev, lumi_ramp_time):
        self.bunches = bunches
        self.interactions = interactions
        self.lumi_start = lumi_start
        self.lumi_lev = lumi_lev
        self.lumi_ramp_time = lumi_ramp_time





class LumiCalculator:
    proton_mass = 0.938 # in GeV
    emitt = 2.5e-6 # [m]
    frev =  299792458/26658.8832# hz
    t1 = 0 #-50e-12# bunch time delay [s]
    t2 = 0 #10e-12# bunch time delay [s]
    #wcc = 400.0 * 1e6 # RF frequency ???


    def __init__(self, crossing_angle, crab_angle, energy, nbunches, nparticles, betx, bety, sig_s, wcc=400*1e6):
        self.wcc = wcc
        self.crab_angle = crab_angle # in radians
        self.crossing_angle = crossing_angle # rad
        self.Nb = nbunches
        self.Np = nparticles
        self.gamma = energy/self.proton_mass
        self.geom_emitt = self.emitt/self.gamma
        
        self.betx = betx
        self.bety = bety
        
        self.sigx = np.sqrt(self.geom_emitt*self.betx)
        self.sigy = np.sqrt(self.geom_emitt*self.bety)
        self.sig_s = sig_s



    


