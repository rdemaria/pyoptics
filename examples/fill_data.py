t1=optics.open('twiss_lhcb1.tfs')
t2=optics.open('twiss_lhcb2.tfs')

basedir='data_ats_md2/betabeat/7-5-2011/'
fnames=sorted(glob(basedir+'LHCB1/Measurements/*/*.gz'))

fn=fnames[-2]
m=betabeatBPM(fn)
m.mk_fitlsq()


m1=betabeatBPM(fnames[-4])
m2=betabeatBPM(fnames[-2])

m1.mk_fitlsq()
m2.mk_fitlsq()

f1=m1.fit_tbl.filter_tune(0.27,0.29)
f2=m2.fit_tbl.filter_tune(0.27,0.29)

plot(f1.s,   f1.co-f2.co)


t=m.fit_tbl
for i in where(t.res<20)[0]:
  figure()
  plot(m.data[t.iii[i]])
  co,ff,a,p=t.co[i],t.tune[i],t.amp[i],t.phase[i]
  plot(ham(m.t,ff,a,p)+co)
  title('%s  %10.5f'%( m.bpms[i], t.res[i]))


figure()
hist(m.fit_tbl.tune,bins=1000)

m.refine_fit()

m.plot_betx(t1)
m.plot_bety(t1)

m.plot_mux(t1)
m.plot_muy(t1)

s=m.fitx.s[1:]
plot(s,-rad2deg(diff(unwrap(m.fitx.phase))))

ix1=array([ where(t1.name==bp)[0][0] for bp in m.fitx.bpms ])
muA=t1.mux[ix1]

plot(s,-rad2deg(diff(unwrap(muA*2*pi))))

iy1=array([ where(t1.name==bp)[0][0] for bp in m.fity.bpms ])

betxA=sqrt(t1.betx[ix1])
betxB=m.fitx.amp
fact=(betxA/betxB).mean()
#fact=optimize.fmin(lambda x: sum((betxA-x[0]*betxB)**2),[fact])

plot(m.fitx.s,betxA,'b-')
plot(m.fitx.s,betxB*fact,'r-')


run __init__.py
m.__class__=betabeatBPM
m.refine_fit()



t=arange(m.turns)
f=t2f(t)


out=[]
for i in range(len(m.data)):
  if not m.bad[i]:
    ff,a,p,res=lsqmax(m.data[i],t)
    s=m.spos[i]
    out.append([s,ff,a,p,res])


outx=[]
outy=[]
for s,ff,a,p,res in out:
    if m.flags[i]==0:
      if abs(ff-0.27)<0.01:
        outx.append([s,ff,a,p,res])
    else:
      if abs(ff-0.32)<0.01:
        outy.append([s,ff,a,p,res])

sss,fff,aaa,ppp=array(sorted(outx)).T

sss,fff,aaa,ppp=array(sorted(outy)).T







for i in range(len(m.data)):
  if not m.bad[i]:
    v=m.data[i]
    fv=rfft(v)
    if any(fv==0):
      print i,m.bpms[i]




for fn in fnames:
  m=betabeatBPM(fn)
  print sum(m.bad)



    plot(f,c2db(rfft(v)))

for v in m.data:
  print abs(diff(v)).min()




v=m.data[0]
plot(t,v)
ff,a,p=getfftmax(v,t2f(t))
plot(t,ham(t,ff,a,p))
ff,a,p=lsqmax(v,t)
plot(t,ham(t,ff,a,p))



print ff,a,p
plot(t,ham(t,ff,a,p))

t2=t[:-1]
v2=ham(t2,0.4,2,0)
ff,a,p=getfftmax(v2,t2f(t2))
plot(t2,v2)
plot(t2,ham(t2,ff,a,p))



for v in m.data[~m.bad]:
  fv=c2db(rfft(v))
  fv=0

  peaks=sorted(zip(fv,range(len(fv))),reverse=True)[:5]


  plot(f,c2db(rfft(v)))




fnames=sorted(glob('data_ats_md2/fill_data/*/BPM/Beam1@Turn@*.sdds.gz'))


fn=fnames[-6]

m=LHCBPM(fn)


run __init__.py
m.__class__=LHCBPM
m.plot_2dtunes()


t1=optics.open('twiss_lhcb1.tfs')
t2=optics.open('twiss_lhcb2.tfs')
m=LHCBPM('data_ats_md2/fill_data/1769/BPM/Beam2@Turn@2011_05_07@16_21_30_415.sdds.gz')

m=LHCBPM('data_ats_md2/fill_data/1769/BPM/Beam1@Turn@2011_05_07@16_27_44_105.sdds.gz')
m=LHCBPM('data_ats_md2/fill_data/1770/BPM/Beam1@Turn@2011_05_07@18_20_06_686.sdds.gz')

m=LHCBPM('data_ats_md2/fill_data/1770/BPM/Beam2@BunchTurn@2011_05_07@18_53_42_682.sdds.gz')

m.mk_spos(t2)
m.mk_fitlsq()
m.mk_fitlsq2()
m.plot_2dtunes()

m=LHCBPM('data_ats_md2/fill_data/1768/BPM//Beam2@Turn@2011_05_07@13_02_14_037.sdds.gz')

