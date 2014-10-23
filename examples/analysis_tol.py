from get_tol import *

fn1='/afs/cern.ch/eng/lhc/optics/V6.503/aperture/aper_tol.b1.madx'
fn2='/afs/cern.ch/eng/lhc/optics/V6.503/aperture/aper_tol.b2.madx'

ttol=read_aptol(fn1)
ttol.update(read_aptol(fn2))

for name in 'MBRC','MQY.4[LR][15]','MQM','MQML':
  summ_tol(ttol,name)


