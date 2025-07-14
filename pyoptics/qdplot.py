#from __future__ import print_function

import re
import os
import gzip
import time


import matplotlib.pyplot as pl
import matplotlib
from matplotlib.animation import FuncAnimation
import numpy as np
import scipy


#from .utils import mystr as _mystr
#from .utils import pyname
#from collections import namedtuple
#from .pydataobj import dataobj
#from . import tfsdata
#from .survey import rot_mad, get_w_from_angles
#from .tablemixin import TableMixIn


def _mylbl(d, x):
    return d.get(x, r"$%s$" % x)


class qdplot(object):
    lglabel = {
        "betx": r"$\beta_x$",
        "bety": r"$\beta_y$",
        "dx": r"$D_x$",
        "dy": r"$D_y$",
        "mux": r"$\mu_x$",
        "muy": r"$\mu_y$",
        "Ax": "$A_x$",
        "Ay": "$A_y$",
        "Bx": "$B_x$",
        "By": "$B_y$",
        "wx": "$w_x$",
        "wy": "$w_y$",
        "sigx": r"$\sigma_x=\sqrt{\beta_x \epsilon}$",
        "sigy": r"$\sigma_y=\sqrt{\beta_y \epsilon}$",
        "sigdx": r"$\sigma_{D_x}=D_x \delta$",
        "n1": r"Aperture [$\sigma$]",
    }

    axlabel = {
        "s": r"$s [m]$",
        "ss": r"$s [m]$",
        "betx": r"$\beta [m]$",
        "bety": r"$\beta [m]$",
        "mux": r"$\mu/(2 \pi)$",
        "muy": r"$\mu/(2 \pi)$",
        "dx": r"$D [m]$",
        "dy": r"$D [m]$",
        "x": r"$co [m]$",
        "y": r"$co [m]$",
        "sigx": r"$\sigma$ [mm]",
        "sigy": r"$\sigma$ [mm]",
        "sigdx": r"$\sigma$ [mm]",
        "n1": r"Aperture [$\sigma$]",
    }
    autoupdate = []

    def ani_autoupdate(self):
        from matplotlib.animation import FuncAnimation

        self._ani = FuncAnimation(self.figure, self.update, blit=False, interval=1000)

    def ani_stopupdate(self):
        del self._ani

    @classmethod
    def on_updated(cls, fun):
        cls.on_update = fun

    def __init__(
        self,
        t,
        x="",
        yl="",
        yr="",
        idx=slice(None),
        clist="C1 C2 C3 C4 C5 C6 C7 C8",
        lattice=None,
        newfig=True,
        pre=None,
    ):
        yl, yr, clist = list(map(str.split, (yl, yr, clist)))
        #    timeit('Init',True)
        self.color = {}
        self.left = None
        self.right = None
        self.lattice = None
        self.pre = None
        self.t, self.x, self.yl, self.yr, self.idx, self.clist = (
            t,
            x,
            yl,
            yr,
            idx,
            clist,
        )
        for i in self.yl + self.yr:
            self.color[i] = self.clist.pop(0)
            self.clist.append(self.color[i])
        if newfig is True:
            self.figure = pl.figure()
        elif newfig is False:
            self.figure = pl.gcf()
            self.figure.clf()
        else:
            self.figure = newfig
            self.figure.clf()
        if lattice:
            self.lattice = self._new_axes()
            #      self.lattice.set_autoscale_on(False)
            self.lattice.yaxis.set_visible(False)
        if yl:
            self.left = self._new_axes()
            #      self.left.set_autoscale_on(False)
        if yr:
            self.right = self._new_axes()
            #      self.right.set_autoscale_on(False)
            self.left.yaxis.set_label_position("right")
            self.left.yaxis.set_ticks_position("right")

        #    timeit('Setup')
        self.run()
        if lattice is not None:
            self.lattice.set_autoscale_on(False)
        if yl:
            self.left.set_autoscale_on(False)
            self.left.yaxis.set_label_position("left")
            self.left.yaxis.set_ticks_position("left")
        if yr:
            self.right.set_autoscale_on(False)
            self.right.yaxis.set_label_position("right")
            self.right.yaxis.set_ticks_position("right")

    #    timeit('Update')
    def _new_axes(self):
        if self.figure.axes:
            ax = self.figure.axes[-1]
            out = self.figure.add_axes(ax.get_position(), sharex=ax, frameon=False)
        else:
            # adjust plot dimensions
            out = self.figure.add_axes([0.17, 0.12, 0.6, 0.8])
        return out

    def __repr__(self):
        return object.__repr__(self)

    def _trig(self):
        print("optics trig")
        self.run()

    def update(self, *args):
        if hasattr(self.t, "reload"):
            if self.t.reload():
                self.run()
                return self
        return False

    #  def _wx_callback(self,*args):
    #    self.update()
    #    wx.WakeUpIdle()
    #
    #  def autoupdate(self):
    #    if pl.rcParams['backend']=='WXAgg':
    #      wx.EVT_IDLE.Bind(wx.GetApp(),wx.ID_ANY,wx.ID_ANY,self._wx_callback)
    #    return self
    #
    #  def stop_update(self):
    #    if pl.rcParams['backend']=='WXAgg':
    #      wx.EVT_IDLE.Unbind(wx.GetApp(),wx.ID_ANY,wx.ID_ANY,self._callback)
    #
    #  def __del__(self):
    #    if hasattr(self,'_callback'):
    #      self.stop_update()

    def run(self):
        #    print 'optics run'
        self.ont = self.t
        self.xaxis = getattr(self.ont, self.x)[self.idx]
        is_ion = pl.isinteractive()
        pl.interactive(False)
        self.lines = []
        self.legends = []
        #    self.figure.lines=[]
        #    self.figure.patches=[]
        #    self.figure.texts=[]
        #    self.figure.images = []
        self.figure.legends = []

        if self.lattice:
            self.lattice.clear()
            self._lattice(["k0l", "kn0l", "angle"], "#a0ffa0", "Bend h")
            self._lattice(["ks0l"], "#ffa0a0", "Bend v")
            self._lattice(["kn1l", "k1l"], "#a0a0ff", "Quad")
            self._lattice(["hkick"], "#e0a0e0", "Kick h")
            self._lattice(["vkick"], "#a0e0e0", "Kick v")
            self._lattice(["kn2l", "k2l"], "#e0e0a0", "Sext")
        if self.left:
            self.left.clear()
            for i in self.yl:
                self._column(i, self.left, self.color[i])
        if self.right:
            self.right.clear()
            for i in self.yr:
                self._column(i, self.right, self.color[i])
        ca = self.figure.gca()
        ca.set_xlabel(_mylbl(self.axlabel, self.x))
        ca.set_xlim(min(self.xaxis), max(self.xaxis))
        self.figure.legend(self.lines, self.legends, loc="upper right")
        ca.grid(True)
        #    self.figure.canvas.mpl_connect('button_release_event',self.button_press)
        self.figure.canvas.mpl_connect("pick_event", self.pick)
        pl.interactive(is_ion)
        self.figure.canvas.draw()
        if hasattr(self, "on_run"):
            self.on_run(self)

    def pick(self, event):
        pos = np.array([event.mouseevent.x, event.mouseevent.y])
        name = event.artist.elemname
        prop = event.artist.elemprop
        value = event.artist.elemvalue
        print("\n %s.%s=%s" % (name, prop, value), end=" ")

    #  def button_press(self,mouseevent):
    #    rel=np.array([mouseevent.x,mouseevent.y])
    #    dx,dy=self.pickpos/rel
    #    print 'release'
    #    self.t[self.pickname][self.pickprop]*=dy
    #    self.t.track()
    #    self.update()

    def _lattice(self, names, color, lbl):
        #    timeit('start lattice %s' % names,1)
        vd = 0
        sp = self.lattice
        s = self.ont.s
        l = self.ont.l
        for i in names:
            myvd = self.ont.__dict__.get(i, None)
            if myvd is not None:
                vdname = i
                vd = myvd[self.idx] + vd
        if np.any(vd != 0):
            m = np.abs(vd).max()
            if m > 1e-10:
                c = np.where(abs(vd) > m * 1e-4)[0]
                if len(c) > 0:
                    if np.all(l[c] > 0):
                        vd[c] = vd[c] / l[c]
                        m = abs(vd[c]).max()
                    vd[c] /= m
                    if self.ont._is_s_begin:
                        plt = self.lattice.bar(
                            s[c] + l[c] / 2, vd[c], l[c], picker=True
                        )  # changed
                    else:
                        plt = self.lattice.bar(
                            s[c] - l[c] / 2, vd[c], l[c], picker=True
                        )  # changed
                    pl.setp(plt, facecolor=color, edgecolor=color)
                    if plt:
                        self.lines.append(plt[0])
                        self.legends.append(lbl)
                    row_names = self.ont.name
                    for r, i in zip(plt, c):
                        r.elemname = row_names[i]
                        r.elemprop = vdname
                        r.elemvalue = getattr(self.ont, vdname)[i]
                self.lattice.set_ylim(-1.5, 1.5)

    #    timeit('end lattice')

    def _column(self, name, sp, color):
        fig, s = self.figure, self.xaxis
        y = self.ont(name)[self.idx]
        (bxp,) = sp.plot(s, y, color, label=_mylbl(self.lglabel, name))
        sp.set_ylabel(_mylbl(self.axlabel, name))
        self.lines.append(bxp)
        self.legends.append(_mylbl(self.lglabel, name))
        sp.autoscale_view()

    def savefig(self, name):
        self.figure.savefig(name)
        return self


