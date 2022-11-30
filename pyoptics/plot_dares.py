from numpy import *
from matplotlib.pyplot import *
from glob import glob
import os
import numpy as np


def dict_to_array(data, names=["id", "data"], formats=[object, object]):
    ftype = dict(names=names, formats=formats)
    return np.array(data.items(), dtype=ftype)


def load_dares(fh):
    out = []
    for l in fh:
        l = l.split()
        out.append(map(float, l[1:]))
    return [array(i) for i in zip(*out)]


def plot_stat(fname, lcolor="m", lbl="", ymin=8, ymax=20):
    data = load_dares(open(fname))
    if len(data) == 7:
        angle, damin, daavg, damax, daflag, daini, daout = data
    else:
        angle, damin, daavg, damax, daflag, daini, daout, dastd = data
    angle *= 90 / float(angle[-1] + 1)
    plot(angle, abs(damin), linestyle="--", color=lcolor)
    plot(angle, abs(daavg), linestyle="-", color=lcolor, label=lbl)
    plot(angle, abs(damax), linestyle="--", color=lcolor)
    grid(True)
    xlabel(r"Angle [degree]")
    ylabel(r"DA $[\sigma]$")
    xlim(0, 90)
    ylim(ymin, ymax)


#  return data


def plot_stat_avg(fname, color="m", lbl=""):
    data = load_dares(open(fname))
    angle, damin, daavg, damax, daflag, daini, daout = data
    angle *= 90 / float(angle[-1] + 1)
    plot(angle, abs(daavg), "-%s" % color, label=lbl)
    xlim(0, 90)
    grid(True)
    xlabel(r"Angle [degree]")
    ylabel(r"Avg DA $[\sigma]$")
    ylim(8, 20)


#  return data


def plot_stat_min(fname, color="m", lbl=""):
    data = load_dares(open(fname))
    angle, damin, daavg, damax, daflag, daini, daout = data
    angle *= 90 / float(angle[-1] + 1)
    plot(angle, abs(damin), "-%s" % color, label=lbl)
    xlim(0, 90)
    grid(True)
    xlabel(r"Angle [degree]")
    ylabel(r"Min DA $[\sigma]$")
    ylim(8, 20)


#  return data


def plot_stat_std(fname, lcolor="m", lbl="", ymin=8, ymax=20, nsigma=1):
    data = load_dares(open(fname))
    angle, damin, daavg, damax, daflag, daini, daout, dastd = data
    angle *= 90 / float(angle[-1] + 1)
    plot(angle, abs(damin), linestyle="--", color=lcolor)
    plot(angle, abs(daavg), linestyle="-", color=lcolor, label=lbl)
    plot(angle, abs(damax), linestyle="--", color=lcolor)
    fill_between(
        angle,
        daavg - nsigma * dastd,
        daavg + nsigma * dastd,
        facecolor=lcolor,
        alpha=0.25,
    )
    grid(True)
    xlabel(r"Angle [degree]")
    ylabel(r"DA $[\sigma]$")
    xlim(0, 90)
    ylim(ymin, ymax)


#  return data


def plot_all_seeds(
    dname, lcolor="m", lbl="", ymin=8, ymax=20, dacrit="daavgamp", lseeds=None
):
    """plots the DA vs angle for each seed,
    dname is the directory name of the DAres files dares_*
    lseeds=list of seeds e.g. [1,4,5,20]
    daslope = Strict chaotic boundary via slope method
    daspace = Certain chaotic boundary via large distance in phase space method
    daavgamp= Dynamic aperture concerning the phase space averaged amplitude
    dainamp = Raw dynamic aperture concerning initial amplitude
    daini   = Lower bound of tracked amplitude range
    daout   = Upper bound of tracked amplitude range"""
    close("all")
    idx = {
        "daslope": 0,
        "daspace": 1,
        "daavgamp": 3,
        "dainamp": 4,
        "daini": 5,
        "daout": 6,
    }
    ida = idx[dacrit]
    data = {}
    for fname in glob(dname + "/DAres*[0-9]"):
        angle = fname.split(".")[-1]
        data[angle] = load_dares(open(fname))
    nseeds = len(data["1"][0])  # number of seeds
    nangles = len(data.keys())
    dangle = dict.fromkeys(range(nseeds), [])
    dda = dict.fromkeys(range(nseeds), [])
    for seed in dangle.keys():
        for angle in sorted(data.keys(), key=int):
            dangle[seed] = dangle[seed] + [float(angle) * 90 / (nangles + 1)]
            dda[seed] = dda[seed] + [data[angle][ida][seed]]
    if lseeds == None:
        lseeds = dangle.keys()
    for seed in lseeds:  # list of seeds starts with 1, while index starts with 0
        plot(
            dangle[seed - 1], abs(np.array(dda[seed - 1])), linestyle="--", color=lcolor
        )
    grid(True)
    xlabel(r"Angle [degree]")
    ylabel(r"DA $[\sigma]$")
    xlim(0, 90)
    ylim(ymin, ymax)
