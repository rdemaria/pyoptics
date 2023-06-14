from pyoptics import lumi2

#HL-LHC
ip1=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetay=250e-6*2,ccy=380e-6)
ip2=lumi2.InteractionPoint(betx=10,bety=10,thetay=200e-6*2)
ip5=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetax=250e-6*2)
ip8=lumi2.InteractionPoint(betx=1.5,bety=1.5,thetay=180e-6*2)
bunch=lumi2.Bunch(ips=(ip1,ip5,ip2,ip8))

print(ip1.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
print(ip1.lumi(bunch)/1e38, "cm^-2 s^-1")
print(ip1.instantaneous_lumi(bunch)/1e38, "cm^-2 s^-1")

print(ip2.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
print(ip2.lumi(bunch)/1e38, "cm^-2 s^-1")
print(ip2.instantaneous_lumi(bunch)/1e38, "cm^-2 s^-1")

# bunch.get_lumi()

#LHC
b_lhc=lumi2.Bunch(energy=7e12,emitnx=3.75e-6,emitny=3.75e-6,nb=2808,ppb=1.1e11)
ip1_lhc=lumi2.InteractionPoint(betx=0.55,bety=0.55, thetax=285e-6)
ip1_lhc.lumi(b_lhc)/1e38


def test_rf():
    print("---------------With CC RF---------------")
    print("Luminosity for at collision adjustment")
    # ip1=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetay=250e-6*2,ccy=380e-6, ccrf=400e6)
    # ip2=lumi2.InteractionPoint(betx=10,bety=10,thetay=200e-6*2)
    # ip5=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetax=250e-6*2,ccx=380e-6, ccrf=400e6)
    # ip8=lumi2.InteractionPoint(betx=1.5,bety=1.5,thetay=180e-6*2)
    betxy_i5 = 1
    ip1=lumi2.InteractionPoint(betx=betxy_i5,bety=betxy_i5,thetay=250e-6*2,ccy=380e-6, ccrf=400.79e6)
    ip2=lumi2.InteractionPoint(betx=10,bety=10,thetay=200e-6*2)
    ip5=lumi2.InteractionPoint(betx=betxy_i5,bety=betxy_i5,thetax=250e-6*2,ccx=380e-6, ccrf=400.79e6)
    ip8=lumi2.InteractionPoint(betx=1.5,bety=1.5,thetay=180e-6*2)
    bunch=lumi2.Bunch(ips=(ip1,ip5,ip2,ip8))
    print("IP1")
    print(ip1.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
    print(ip1.lumi(bunch)/1e38, "cm^-2 s^-1")
    print(ip1.instantaneous_lumi(bunch)/1e38, "cm^-2 s^-1")
    print("IP5")
    print(ip5.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
    print(ip5.lumi(bunch)/1e38, "cm^-2 s^-1")
    print(ip5.instantaneous_lumi(bunch)/1e38, "cm^-2 s^-1")

    bunch=lumi2.Bunch(ips=(ip2,ip8))
    print("IP2")
    print(ip2.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
    print(ip2.lumi(bunch)/1e38, "cm^-2 s^-1")

    print("IP8")
    print(ip8.lumi_headon(bunch)/1e38, "cm^-2 s^-1")
    print(ip8.lumi(bunch)/1e38, "cm^-2 s^-1")
