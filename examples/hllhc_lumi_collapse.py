from pyoptics import lumi


### collapse
ip1 = lumi.IP(name="atlas", betx=1, bety=1, py=250e-6)
b15 = lumi.Bunch(nb=2748, ppb=2.2e11, ips=(ip1,))

print(b15.luminosity(verbose=False) / 1e38, "10^34 cm^-2 s^-1")

print(f"nb   energy  lumi   pileup  betastar")
for target in [1, 2, 2.5]:
    for energy in [6.8, 7]:
        for nb in [1960, 2200, 2748]:
            b15c = b15.clone(nb=nb, energy=energy * 1e12)
            ipc = ip1.betastar_from_lumi(bunch=b15c, target=target * 1e38)
            ll = ipc.luminosity(bunch=b15c) / 1e38
            pileup = ipc.pileup(bunch=b15c)
            print(
                f"{nb:4}   {energy:4.1f}   {ll:4.1f}     "
                f"{pileup:4.0f}        {ipc.betx:4.2f}"
            )

ipc=ip1.clone(betx=0.5,bety=0.2)
print(f"nb   energy  lumi   pileup  betastar")
for target in [2.5]:
    for energy in [6.8, 7]:
        for nb in [1960, 2200, 2748]:
            b15c = b15.clone(nb=nb, energy=energy * 1e12)
            ipc = ipc.betastar_from_lumi(bunch=b15c, target=target * 1e38, betaratio=4)
            ll = ipc.luminosity(bunch=b15c) / 1e38
            pileup = ipc.pileup(bunch=b15c)
            print(
                f"{nb:4}   {energy:4.1f}   {ll:4.1f}     "
                f"{pileup:4.0f}        {ipc.betx:4.2f}  {ipc.bety:4.2f}"
            )
