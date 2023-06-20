from pyoptics import lumi


ip8 = lumi.IP(
    name="lhcb",
    betx=1.5,
    bety=1.5,
    px=135e-6,
    py=170e-6 + 1.8e-6 )

b8=lumi.Bunch(nb=2572, ips=(8))

ip8.clone(px=135e-6,py=1.8e-6+160e-6).luminosity(b8) / 1e38
ip8.clone(px=135e-6-200e-6,py=1.8e-6).luminosity(b8)/1e38
ip8.clone(px=135e-6+150e-6,py=1.8e-6).luminosity(b8)/1e38


ip8.clone(px=135e-6,py=1.8e-6+170e-6).luminosity(b8) / 1e38
