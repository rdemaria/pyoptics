import sddsdata

from matplotlib.pyplot import *
from numpy import *

from harmonic_fit import *
from rdmdate import *


class LHCBPM(object):
    def __init__(self,fn):
        self.fn=fn
        sdds=sddsdata.sddsdata(open(fn))
        self.bpms=r_[sdds.data[0]['bpmNames']]
        self.nbpms=len(self.bpms)
        self.turns=sdds.data[0]['nbOfCapTurns'][0]
        self.bunches=sdds.data[0]['nbOfCapBunches'][0]
        xdata=sdds.data[0]['horPositionsConcentratedAndSorted']
        ydata=sdds.data[0]['verPositionsConcentratedAndSorted']
        self.xbpm=xdata.reshape(self.nbpms,self.bunches,self.turns)
        self.ybpm=ydata.reshape(self.nbpms,self.bunches,self.turns)
        self.dt=sdds.data[0]['acqStamp']/1e9
        self.sdds=sdds
        #print self.nbpms,self.bunches,self.turns
    def get_xy(self,bpm,bunch):
        return self.xbpm[bpm,bunch,:],self.ybpm[bpm,bunch,:]
    def plot_xy(self,bpm,bunch):
        x,y=self.get_xy(bpm,bunch)
        plot(x,label='x[b=%d] %s'%(bunch,self.bpms[bpm],))
        plot(y,label='y[b=%d] %s'%(bunch,self.bpms[bpm],))
        title(dumpdate(self.dt))
        legend()
    def plot_xy_fft(self,bpm,bunch):
        x,y=self.get_xy(bpm,bunch)
        f=linspace(0,0.5,len(x)/2+1)
        plot(f,abs(rfft(x)),label='x[b=%d] %s'%(bunch,self.bpms[bpm],))
        plot(f,abs(rfft(y)),label='y[b=%d] %s'%(bunch,self.bpms[bpm],))
        title(dumpdate(self.dt))
        legend()
    def get_tune(bpm,bunch,tune_fit=maxharm_brent):
        x,y=self.get_xy(bpm,bunch)
        return tune_fit(x-x.mean())[0],tune_fit(y-y.mean())[0]
    def get_tunes(self,tune_fit=maxharm_brent):
        qx=[];qy=[]
        x=self.xbpm;y=self.ybpm
        bpms,bunches,turns=x.shape
        for b in  range(self.bunches):
          for p in range(self.nbpms):
            xx=x[p,b,:]; yy=y[p,b,:];
            if sum(xx**2)>0:
                qx.append(tune_fit(xx-xx.mean(),0.20,0.40))
            if sum(yy**2)>0:
                qy.append(tune_fit(yy-yy.mean(),0.20,0.40))
        return array(qx),array(qy)
    def plot_tune_hist(self):
        tunex,tuney=self.get_tunes()
        qx,ax,px,rx=zip(*tunex)
        qy,ay,py,ry=zip(*tuney)
        hist(qx,bins=10,label='Qx');
        hist(qy,bins=10,label='Qy');
        legend()
        xlabel('Tune')
        ylabel('Count')
        title(dumpdate(self.dt))
    def get_coupling(self):
        out=[]
        for b in  range(self.bunches):
          for p in range(self.nbpms):
            print b,p
            xx=self.xbpm[p,b,:]; yy=self.ybpm[p,b,:]
            cc=fit_coupled_lsq2(xx,yy)
            out.append([b,p]+list(cc))
        return out


