from pyoptics import lumi


### start of levelling
ip1 = lumi.IP(name="atlas", betx=0.5, bety=0.5, py=250e-6, ccy=-190e-6)
b15 = lumi.Bunch(nb=2748, ppb=2.2e11, ips=(ip1,))

print(b15.luminosity(verbose=True) / 1e38, "10^34 cm^-2 s^-1")

print(f"nb   energy  lumi   pileup  betastar")
for target in [140, 200]:
    for energy in [6.8, 7]:
        b15c = b15.clone(nb=nb, energy=energy * 1e12)
        ipc = ip1.betastar_from_pileup(bunch=b15c, target=target)
        pileup = ipc.pileup(bunch=b15c)
        for nb in [1960, 2200, 2748]:
            b15c = b15.clone(nb=nb, energy=energy * 1e12)
            ll = ipc.luminosity(bunch=b15c) / 1e38
            print(
                f"{nb:4}   {energy:4.1f}   {ll:4.1f}     "
                f"{pileup:4.0f}        {ipc.betx:4.2f}"
            )

print(f"nb   energy  lumi   pileup  betastar")
for target in [140, 200]:
    for energy in [6.8, 7]:
        b15c = b15.clone(nb=nb, energy=energy * 1e12)
        ipc = ip1.betastar_from_pileup(bunch=b15c, target=target, betaratio=4)
        pileup = ipc.pileup(bunch=b15c)
        for nb in [1960, 2200, 2748]:
            b15c = b15.clone(
                nb=nb,
                energy=energy * 1e12,
                betx=0.5 * 140 / target,
                bety=0.5 * 140 / target*4,
            )
            ll = ipc.luminosity(bunch=b15c) / 1e38
            print(
                f"{nb:4}   {energy:4.1f}   {ll:4.1f}     "
                f"{pileup:4.0f}        {ipc.betx:4.2f}  {ipc.bety:4.2f}"
            )

