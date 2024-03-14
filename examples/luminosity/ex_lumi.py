from pyoptics import lumi


# Run4 Start of collapse
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

b15 = lumi.Bunch(nb=2748, ips=(ip1, ip5))
b2 = b15.clone(nb=2492, ips=(ip2,))
b8 = b15.clone(nb=2574, ips=(ip8,))
ip8 = ip8.sep_from_lumi(b8, 2e37)
ip2 = ip2.sep_from_lumi(b2, 1.4e35)


print(ip1.luminosity(b15, verbose=False) / 1e38, "10^34 cm^-2 s^-1")
print(ip5.luminosity(b15, verbose=False) / 1e38, "10^34 cm^-2 s^-1")
print(ip2.luminosity(b2, verbose=False) / 1e38, "10^34 cm^-2 s^-1")
print(ip8.luminosity(b8, verbose=False) / 1e38, "10^34 cm^-2 s^-1")

for nb in [2748, 2200, 1960]:
    for energy in [6.8e12, 7e12]:
        b15c = b15.clone(nb=nb, energy=energy)
        for target in [1e38, 2e38, 2.5e38]:
            ipc = ip1.solve_luminosity_roundbeam(bunch=b15c, target=target)
            print(f"nb={nb} energy={energy} target={target} betastar={ipc.betx}")


# Lumi levelling
ip1 = lumi.IP(name="atlas", betx=0.2, bety=0.2, py=250e-6, ccy=-190e-6)
ip2 = lumi.IP(
    name="alice", betx=10, bety=10, py=-70e-6 - 170e-6, px=0.3e-6, sepx=0.135e-3
)
ip5 = lumi.IP(name="cms", betx=0.2, bety=0.2, px=250e-6, ccx=-190e-6)
ip8 = lumi.IP(
    name="lhcb",
    betx=1.5,
    bety=1.5,
    px=135e-6,
    py=170e-6 + 1.8e-6,
    sepx=0.032e-3,
)
b15 = lumi.Bunch(
    nb=2748, ppb=1.3e11, emitnx=2.9e-6, emitny=2.2e-6, ips=(ip1, ip5, ip2, ip8)
)
b2 = b15.clone(nb=2492, ips=(ip2,))
b8 = b15.clone(nb=2574, ips=(ip8,))

print(b15.luminosity(verbose=False) / 1e38, "10^34 cm^-2 s^-1")
print(b2.luminosity(verbose=False) / 1e38, "10^34 cm^-2 s^-1")
print(b8.luminosity(verbose=False) / 1e38, "10^34 cm^-2 s^-1")


# HL-LHC nominal
ip1 = lumi.IP(name="atlas", betx=0.15, bety=0.15, py=250e-6, ccy=-190e-6)
ip5 = lumi.IP(name="cms", betx=0.15, bety=0.15, px=250e-6, ccx=-190e-6)
b15 = lumi.Bunch(nb=2748, ppb=2.2e11, emitnx=2.5e-6, emitny=2.5e-6, ips=(ip1, ip5))
ip1.luminosity(b15, verbose=True) / 1e38
ip5.luminosity(b15, verbose=True) / 1e38
ip1.luminosity(b15) / 1e38
ip5.luminosity(b15) / 1e38
ip1.lumi_simple(b15) / 1e38
ip5.lumi_simple(b15) / 1e38


# LHC

b_lhc = lumi.Bunch(nb=2808, energy=7e12, emitnx=3.75e-6, emitny=3.75e-6, ppb=1.15e11)
ip_lhc = lumi.IP(name="lhc", betx=0.55, bety=0.55, px=142.5e-6)
print(ip_lhc.luminosity(b_lhc, verbose=True) / 1e38, "10^34 cm^-2 s^-1")


# Run4 table 13
ip1 = lumi2.InteractionPoint(
    betx=1, bety=1, thetay=250e-6 * 2, ccy=190e-6 * 2, ccrf=400.79e6 * 0
)
ip2 = lumi2.InteractionPoint(
    betx=10, bety=10, thetay=-100e-6 * 2, thetax=0.3e-6 * 2, sepx=0.13e-3 * 2
)
ip5 = lumi2.InteractionPoint(
    betx=1, bety=1, thetax=250e-6 * 2, ccx=190e-6 * 2, ccrf=400.79e6 * 0
)
ip8 = lumi2.InteractionPoint(
    betx=1.5,
    bety=1.5,
    thetax=135e-6 * 2,
    thetay=170e-6 * 2 + 1.8e-6 * 2,
    ccx=190e-6 * 2,
    ccy=190e-6 * 2,
)
b15 = lumi2.Bunch(nb=2748, ips=(ip1, ip5))
bunch2 = lumi2.Bunch(nb=2492, ips=(ip2,))
bunch8 = lumi2.Bunch(nb=2574, ips=(ip8,))


print(b15.lumi() / 1e38, "10^34 cm^-2 s^-1")
print(bunch2.lumi() / 1e38, "10^34 cm^-2 s^-1")
print(bunch8.lumi() / 1e38, "10^34 cm^-2 s^-1")


# HL-LHC
ip1 = lumi2.InteractionPoint(betx=0.15, bety=0.15, thetay=250e-6 * 2, ccy=380e-6)
ip2 = lumi2.InteractionPoint(betx=10, bety=10, thetay=200e-6 * 2)
ip5 = lumi2.InteractionPoint(betx=0.15, bety=0.15, thetax=250e-6 * 2)
ip8 = lumi2.InteractionPoint(betx=1.5, bety=1.5, thetay=180e-6 * 2)
bunch = lumi2.Bunch(ips=(ip1, ip5, ip2, ip8))

print(ip1.lumi(bunch) / 1e38, "10^34 cm^-2 s^-1")
print(ip5.lumi(bunch) / 1e38, "10^34 cm^-2 s^-1")
print(ip2.lumi(bunch) / 1e38, "10^34 cm^-2 s^-1")
print(ip8.lumi(bunch) / 1e38, "10^34 cm^-2 s^-1")

# bunch.get_lumi()

# LHC
b_lhc = lumi2.Bunch(energy=7e12, emitnx=3.75e-6, emitny=3.75e-6, nb=2808, ppb=1.1e11)
ip1_lhc = lumi2.InteractionPoint(betx=0.55, bety=0.55, thetax=285e-6)
ip1_lhc.lumi(b_lhc) / 1e38


def test_rf():
    print("---------------With CC RF---------------")
    print("Luminosity for at collision adjustment")
    # ip1=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetay=250e-6*2,ccy=380e-6, ccrf=400e6)
    # ip2=lumi2.InteractionPoint(betx=10,bety=10,thetay=200e-6*2)
    # ip5=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetax=250e-6*2,ccx=380e-6, ccrf=400e6)
    # ip8=lumi2.InteractionPoint(betx=1.5,bety=1.5,thetay=180e-6*2)

    bunch = lumi2.Bunch(ips=(ip1, ip5, ip2, ip8))
    print("IP1")
    print(ip1.lumi_headon(bunch) / 1e38, "cm^-2 s^-1")
    print(ip1.lumi(bunch) / 1e38, "cm^-2 s^-1")
    print(ip1.instantaneous_lumi(bunch) / 1e38, "cm^-2 s^-1")
    print("IP5")
    print(ip5.lumi_headon(bunch) / 1e38, "cm^-2 s^-1")
    print(ip5.lumi(bunch) / 1e38, "cm^-2 s^-1")
    print(ip5.instantaneous_lumi(bunch) / 1e38, "cm^-2 s^-1")

    bunch = lumi2.Bunch(ips=(ip2, ip8))
    print("IP2")
    print(ip2.lumi_headon(bunch) / 1e38, "cm^-2 s^-1")
    print(ip2.lumi(bunch) / 1e38, "cm^-2 s^-1")

    print("IP8")
    print(ip8.lumi_headon(bunch) / 1e38, "cm^-2 s^-1")
    print(ip8.lumi(bunch) / 1e38, "cm^-2 s^-1")
