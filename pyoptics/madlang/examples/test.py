import madlang
reload(madlang)

mad=madlang.open("examples/sps2.seq")
mad=madlang.open("examples/SPS_Q20_thin.seq")
mad.kmba*=2
mad.mba_40370__3.knl[0]*10==mad.kmba
mad.kmba/=2
out,rest=mad.sps.expand_struct()

from sixtracklib.pysixtrack import *

convert={'drift':Drift, 'mult':Multipole, 'cav': Cavity}

out=[ (name,ccc,convert[el.__class__.__name__](*el)) for name,ccc,el in out]




fn='/home/rdemaria/work/hllhc/HLLHCV1.0/opt_inj.madx'

fn='/afs/cern.ch/eng/lhc/optics/HLLHCV1.0/toolkit/macro.madx'
fn='SPS_Q20_thin.seq'
fn='lhc_as-built_db.seq'

a=madlang.open('examples/lhc/lhc_as-built_db.seq')


print a.mqxa_1l1
madlang.open('examples/lhc/opt_inj.madx',a)
print a.mqxa_1l1

sps=madlang.open("examples/sps2.seq")

sps=madlang.open("examples/SPS_Q20_thin.seq")
out,rest=sps.sps.expand_struct()

s=0
for n,el in o[:10]:
   s+=el.l
   if 'at' in el:
     print n,s,el.at


import madlang
sps=madlang.open("examples/SPS_Q20_thin.seq")
out=sps.build_dep()
sps.print_dep('kqf2')


lhc=madlang.open('examples/lhc_as-built_db.seq')
madlang.open('examples/opt_inj.madx',lhc)
lhc.print_dep('on_sep1')

out=lhc.build_dep()





