from numpy import *
from matplotlib.pyplot import *
from glob import glob
import os

def load_dares(fh):
  out=[]
  for l in fh:
    l=l.split()
    out.append( map(float,l[1:]))
  return [array(i) for i in zip(*out)]

def plot_stat(fname,lcolor='m',lbl='',ymin=8,ymax=20):
  data=load_dares(open(fname))
  if(len(data)==7):
    angle,damin,daavg,damax,daflag,daini,daout=data
  else:
    angle,damin,daavg,damax,daflag,daini,daout,dastd=data
  angle*=90/float(angle[-1]+1)
  plot(angle,abs(damin),linestyle='--',color=lcolor)
  plot(angle,abs(daavg),linestyle='-',color=lcolor,label=lbl)
  plot(angle,abs(damax),linestyle='--',color=lcolor)
  grid(True)
  xlabel(r'Angle [degree]')
  ylabel(r'DA $[\sigma]$')
  xlim(0,90)
  ylim(ymin,ymax)
#  return data

def plot_stat_avg(fname,color='m',lbl=''):
  data=load_dares(open(fname))
  angle,damin,daavg,damax,daflag,daini,daout=data
  angle*=90/float(angle[-1]+1)
  plot(angle,abs(daavg),'-%s'%color,label=lbl)
  xlim(0,90)
  grid(True)
  xlabel(r'Angle [degree]')
  ylabel(r'Avg DA $[\sigma]$')
  ylim(8,20)
#  return data

def plot_stat_min(fname,color='m',lbl=''):
  data=load_dares(open(fname))
  angle,damin,daavg,damax,daflag,daini,daout=data
  angle*=90/float(angle[-1]+1)
  plot(angle,abs(damin),'-%s'%color,label=lbl)
  xlim(0,90)
  grid(True)
  xlabel(r'Angle [degree]')
  ylabel(r'Min DA $[\sigma]$')
  ylim(8,20)
#  return data

def plot_stat_std(fname,lcolor='m',lbl='',ymin=8,ymax=20,nsigma=1):
  data=load_dares(open(fname))
  angle,damin,daavg,damax,daflag,daini,daout,dastd=data
  angle*=90/float(angle[-1]+1)
  plot(angle,abs(damin),linestyle='--',color=lcolor)
  plot(angle,abs(daavg),linestyle='-',color=lcolor,label=lbl)
  plot(angle,abs(damax),linestyle='--',color=lcolor)
  fill_between(angle,daavg-nsigma*dastd,daavg+nsigma*dastd,facecolor=lcolor,alpha=0.25)
  grid(True)
  xlabel(r'Angle [degree]')
  ylabel(r'DA $[\sigma]$')
  xlim(0,90)
  ylim(ymin,ymax)
#  return data
