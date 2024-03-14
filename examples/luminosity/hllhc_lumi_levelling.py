from pyoptics import lumi


ip1 = lumi.IP(name="atlas", betx=1, bety=1, py=250e-6)
ip2 = lumi.IP(name="alice", betx=10, bety=10, py=-70e-6 - 170e-6, px=0.3e-6).sep_(
    0.032e-3
)
ip5 = lumi.IP(name="cms", betx=1, bety=1, px=250e-6)
ip8 = lumi.IP(
    name="lhcb",
    betx=1.5,
    bety=1.5,
    px=135e-6,
    py=170e-6 + 1.8e-6,
).sep_(3.6e-05)


b1528 = lumi.Bunch(nb=2494, ips=(1, 5, 2, 8))
b158 = lumi.Bunch(nb=2572 - 2494, ips=(1, 5, 8))
b15 = lumi.Bunch(nb=2748 - 2572, ips=(1, 5))


b2 = b15.clone(nb=2492, ips=(ip1, ip2, ip5))
b8 = b15.clone(nb=2574, ips=(ip8,))
ip8 = ip8.sep_from_lumi(b8, 2e37)
ip2 = ip2.sep_from_lumi(b2, 1.4e35)


pr = lumi.BetaStarLeveling(
    bunches=(b1528, b158, b15),
    ips={1: ip1, 2: ip2, 5: ip5, 8: ip8},
    lumistart=2.5e38,
    lumilev=2e38,
    lumi_ramp_time=600,
    betastar_end=0.15,
    lumi2=0.014e38,
    lumi8=0.2e38,
)
