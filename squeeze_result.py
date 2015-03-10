import os

import matplotlib.pyplot as pl

import tfsdata
from pydataobj import dataobj


class SqueezeRes(dataobj):
    labels={'betas':r'$\beta^*$ [m]',
            'betx_b1err': r'B1: $(\beta_x-\bar\beta_x)/\bar\beta_x$',
            'bety_b1err': r'B1: $(\beta_y-\bar\beta_y)/\bar\beta_y$',
            'betx_b2err': r'B2: $(\beta_x-\bar\beta_x)/\bar\beta_x$',
            'bety_b2err': r'B2: $(\beta_y-\bar\beta_y)/\bar\beta_y$',
            'mux_b1err': r'B1: $\mu_x-\bar\mu_x$',
            'muy_b1err': r'B1: $\mu_y-\bar\mu_y$',
            'mux_b2err': r'B2: $\mu_x-\bar\mu_x$',
            'muy_b2err': r'B2: $\mu_y-\bar\mu_y$',
            'dx_b1err': r'B1: $D_x-\barD_x$',
            'dy_b1err': r'B1: $D_y-\barD_y$',
            'dx_b2err': r'B2: $D_x-\barD_x$',
            'dy_b2err': r'B2: $D_y-\barD_y$',
            'bet':'Beta beating',
            'd':'Dispersion error',
            'mu':'Local phase error',
            'q':'Tune error',
            'dq':'Chromaticity error'
           }
    @classmethod
    def open(cls,fn):
        data=tfsdata.open(fn)
        data['filename']=fn
        return cls(**data)
    def __init__(self,**data):
        self.__dict__.update(data)
        t=self.ttt/100
        self.step=self.aaa*(1-t)+self.bbb*t
    def plot_curve(self,name,x='betas'):
        lbl=self.labels.get(name,name)
        pl.plot(self[x],self[name],label=lbl)
        return self
    def plot_curves(self,name='bet',xy='xy',x='betas'):
        for beam in '12':
            for yy in xy:
                nnn='%s%s_b%serr'%(name,yy,beam)
                self.plot_curve(nnn,x=x)
        lbl=self.labels.get(x,x)
        pl.xlabel(lbl)
        return self
    def plot_all(self,title='',x='betas'):
        fn=os.path.splitext(self.filename)[0]
        for name in 'bet mu d'.split():
            pl.figure()
            lbl=self.labels.get(name,name)
            pl.title(title+lbl)
            self.plot_curves(name,x=x)
            pl.legend(loc=0)
            pl.savefig('%s_%s.png'%(fn,name))
        for name in 'q dq'.split():
            pl.figure()
            lbl=self.labels.get(name,name)
            pl.title(title+lbl)
            self.plot_curves(name,x=x,xy='12')
            pl.legend(loc=0)
            pl.savefig('%s_%s.png'%(fn,name))
        return self




