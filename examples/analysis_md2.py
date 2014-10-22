
t_start='2011-05-07 10:00:00.000'
t_end='2011-05-07 20:00:00.000'

data=cernlogdb.dbget('RPTE.UA23.RB.A12:I_MEAS',t1=t_start,t2=t_end,conf='ldb.conf')


data1=cernlogdb.dbget(bctdc,t1=t_start,t2=t_end,conf='ldb.conf')
data2=cernlogdb.dbget(bctfr,t1=t_start,t2=t_end,conf='ldb.conf')
data3=cernlogdb.dbget(rffreq,t1=t_start,t2=t_end,conf='ldb.conf')


bpmb1_fns=sorted(glob('data_ats_md2/fill_data/*/BPM/Beam1@Turn@*.sdds.gz'))
bpmb2_fns=sorted(glob('data_ats_md2/fill_data/*/BPM/Beam2@Turn@*.sdds.gz'))
bpmb1_t= cernlogdb.t2num([ bpmdata_date(fn) for fn in bpmb1_fns])
bpmb2_t= cernlogdb.t2num([ bpmdata_date(fn) for fn in bpmb2_fns])

vn=data['datavars'][0]
t,v=data[0]
t=cernlogdb.t2num(t)
plot_date(t,v,'k',label=vn)
ylim(0,8000)
twinx()
cernlogdb.plot_data(data2)
ylim(0,2e10)

twinx()
cernlogdb.plot_data(data3)
ylim(0,2e10)

ax=gca()

[  axvline(t,color='c') for t in bpmb1_t ]
[  axvline(t,color='y') for t in bpmb2_t ]



run harmonic_fit.py
run __init__.py
m.__class__=LHCBPM

bpmb1_fns=sorted(glob('data_ats_md2/fill_data/*/BPM/Beam1@Turn@*.sdds.gz'))
bpmb2_fns=sorted(glob('data_ats_md2/fill_data/*/BPM/Beam2@Turn@*.sdds.gz'))
m=LHCBPM(bpmb1_fns[-7])
m.mk_fitlsq()

t1=optics.open('twiss_lhcb1.tfs')
t2=optics.open('twiss_lhcb2.tfs')
m.mk_spos(t1)


goodx=(~m.badxy) & m.xidx
goody=(~m.badxy) & m.yidx

sum((m.tune*m.res)[goodx])/sum(m.res[goodx])
sum((m.tune*m.res)[goody])/sum(m.res[goody])


u,s,v=svd(m.data[~m.badxy])

u,s,v=svd(m.data[goodx])
f=linspace(0,0.5,2250/2+1)
figure();plot(f,abs(rfft(v[:,:5],axis=0)))







