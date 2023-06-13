from pyoptics import lumi2

#HL-LHC
ip1=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetay=250e-6*2,ccy=380e-6)
ip2=lumi2.InteractionPoint(betx=10,bety=10,thetay=200e-6*2)
ip5=lumi2.InteractionPoint(betx=0.15,bety=0.15,thetax=250e-6*2)
ip8=lumi2.InteractionPoint(betx=1.5,bety=1.5,thetay=180e-6*2)
bunch=lumi2.Bunch(ips=(ip1,ip5,ip2,ip8))

ip1.lumi(bunch)/1e38


bunch.get_lumi()

#LHC
b_lhc=lumi2.Bunch(energy=7e12,emitnx=3.75e-6,emitny=3.75e-6,nb=2808,ppb=1.1e11)
ip1_lhc=lumi2.InteractionPoint(betx=0.55,bety=0.55, thetax=285e-6)
ip1_lhc.lumi(b_lhc)/1e38


