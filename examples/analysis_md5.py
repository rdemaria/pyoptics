from workpkg import *

store=[]

def addtostore(m,name):
  m.name=name
  store.append(m)

fn1=sorted(glob('data_ats_md5/fill_data/FILL_DIR/BPM/Beam1@Turn@*.sdds'))
fn2=sorted(glob('data_ats_md5/fill_data/FILL_DIR/BPM/Beam2@Turn@*.sdds'))
[ '%3d %s'%(i,n.split('/')[-1]) for i,n in enumerate(fn1)]
[ '%3d %s'%(i,n.split('/')[-1]) for i,n in enumerate(fn2)]


m=LHCBPM(fn2[11])
m.mk_spos(optics.open('ats_model10/temp/twiss_lhcb2_12.tfs'))
m.mk_fitlsq()
addtostore(m,'beam 2 1m')

m=LHCBPM(fn1[20])
m.mk_spos(optics.open('ats_model10/temp/twiss_lhcb1_12.tfs'))
m.mk_fitlsq()
addtostore(m,'beam 1 1m')


m=LHCBPM(fn1[27])
m.mk_spos(optics.open('ats_model10/temp/twiss_lhcb1_16.tfs'))
m.mk_fitlsq()
addtostore(m,'beam 1 40cm')
m=LHCBPM(fn2[18])
m.mk_spos(optics.open('ats_model10/temp/twiss_lhcb2_16.tfs'))
m.mk_fitlsq()
addtostore(m,'beam 2 40cm')


m=LHCBPM(fn1[32])
m.mk_spos(optics.open('ats_model10/temp/twiss_lhcb1_16.tfs'))
m.mk_fitlsq()
addtostore(m,'beam 1 40cm corr')



m.plot_beta()


def smallph(m,plane,bpms,bpmt):
  mu=getattr(m,'mu'+plane)
  idx=getattr(m,plane+'idx')
  model=(mu[idx&(m.bpms==bpmt)]-mu[idx&(m.bpms==bpms)])[0]%1
  meas=(m.phase[idx&(m.bpms==bpmt)]-m.phase[idx&(m.bpms==bpms)])[0]%1
  print '%s %s %s mdl:%g meas:%g'%(plane,bpms,bpmt,model,meas)

smallph(m,'x','BPMSW.1L1.B1','BPM.14L1.B1')
smallph(m,'x','BPMSW.1R1.B1','BPM.15R1.B1')
smallph(m,'y','BPMSW.1L1.B1','BPM.15L1.B1')
smallph(m,'y','BPMSW.1R1.B1','BPM.14R1.B1')


smallph(store[2],'y','BPMSW.1L1.B2','BPM.14L1.B2')
smallph(store[2],'y','BPMSW.1R1.B2','BPM.15R1.B2')
smallph(store[2],'x','BPMSW.1L1.B2','BPM.15L1.B2')
smallph(store[2],'x','BPMSW.1R1.B2','BPM.14R1.B2')




def cguess(m1,c,m2):
  phase=arange(0,1,0.001)
  vec=abs(m2-m1*exp(2j*pi*phase)+c)**2
  gphase=phase[vec.argmin()]
  print m2*exp(2j*pi*gphase)

cguess(0.01,0.01j,0.002)

for mm in store:
  mm.__class__=LHCBPM
