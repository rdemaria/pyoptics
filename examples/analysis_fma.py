import os
from glob import glob

from read_fortbin import *

from harmonic_fit import *

#cp /afs/cern.ch/work/m/mfittere/scratch0/w8/track/job_fma_bb_2_2e4/


basedir="/afs/cern.ch/work/m/mfittere/scratch0/w8/track/job_fma_bb_2_2e4/"
sel="1/simul/62.31_60.32/4_6/e4/48/"
fn=glob(os.path.join(basedir,sel,'*.tar'))

head,part=read_allforttar(fn[0])

pdist,x,xp,y,yp,sigma,delta,energy=part[1].T

nfft=8192
out=[]
for i in range(0,len(x)-8192,100):
    print i
    f,a,p,res=maxharm_lsq(y[i:nfft+i])
    out.append([f,a,p,res])

f,a,p,res=zip(*out)
plot(a)

twinx()
plot(delta)

plot(linspace(0,0.5,len(delta)/2+1),log10(abs(rfft(delta))))






