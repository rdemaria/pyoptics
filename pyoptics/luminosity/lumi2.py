class GaussianLongDist:
    def __init__(self,sigz):
        self.sigz = sigz

class QGaussianLongDist:
    def __init__(self,sigz,q):
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
    def __init__(self,name,betx,bety,sepx,sepy,thetax,thetay,dx,dy,ccx,ccy,ccrf,cclag,visible_cross_section=81,total_cross_section=81):
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


class Bunch:
    """Bunch class

    nb: number of particles per bunch
    ppb: number of protons per bunch
    emitx [m.rad]: horizontal emittance
    emity [m.rad]: vertical emittance
    sigz [m]: bunch length
    sigdpp: energy spread
    energy [eV]: beam energy
    ips: list of InteractionPoint objects
    longdist: longitudinal distribution  
    """
    def __init__(self,nb=1976,ppb=2.3e11,emitx=2.5e-6,emity=2.5e-6,sigz=7.5e-3,sigdpp=1.2e-4,energy=7e12,ips=(),longdist=GaussianLongDist):
        self.nb = nb
        self.ppb = ppb
        self.emitx = emitx
        self.emity = emity
        self.sigz = sigz
        self.sigdpp = sigdpp
        self.energy = energy
        self.ips = ips
        self.longdist = longdist


class LevellingProcess:
    def __init__(self, bunches, interactions, lumi_start, lumi_lev, lumi_ramp_time):
        self.bunches = bunches
        self.interactions = interactions
        self.lumi_start = lumi_start
        self.lumi_lev = lumi_lev
        self.lumi_ramp_time = lumi_ramp_time
