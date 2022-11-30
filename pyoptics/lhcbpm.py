from . import sddsdata

from matplotlib.pyplot import *
from numpy import *

from .harmonic_fit_2 import *
from .rdmdate import *


class LHCBPM(object):
    def __init__(self, fn, iqc=True):
        self.fn = fn
        sdds = sddsdata.sddsdata(open(fn))
        self.bpms = r_[sdds.data[0]["bpmNames"]]
        self.nbpms = len(self.bpms)
        self.turns = sdds.data[0]["nbOfCapTurns"][0]
        self.bunches = sdds.data[0]["nbOfCapBunches"][0]
        xdata = sdds.data[0]["horPositionsConcentratedAndSorted"]
        ydata = sdds.data[0]["verPositionsConcentratedAndSorted"]
        if iqc:
            self.xbpm = xdata.reshape(self.nbpms, self.turns, self.bunches).transpose(
                0, 2, 1
            )
            self.ybpm = ydata.reshape(self.nbpms, self.turns, self.bunches).transpose(
                0, 2, 1
            )
        else:
            self.xbpm = xdata.reshape(self.nbpms, self.bunches, self.turns)
            self.ybpm = ydata.reshape(self.nbpms, self.bunches, self.turns)
        self.dt = sdds.data[0]["acqStamp"][0] / 1e9
        self.date = dumpdate(self.dt)
        self.sdds = sdds
        if "B2" in self.bpms[0]:
            self.beam = 2
        elif "B1" in self.bpms[0]:
            self.beam = 1
        # print self.nbpms,self.bunches,self.turns

    def __repr__(self):
        return "<LHCBPM Beam%d %db %s>" % (self.beam, self.bunches, self.date)

    def get_xy(self, bpm, bunch):
        return self.xbpm[bpm, bunch, :], self.ybpm[bpm, bunch, :]

    def plot_xy(self, bpm, bunch):
        x, y = self.get_xy(bpm, bunch)
        plot(
            x,
            label="x[b=%d] %s"
            % (
                bunch,
                self.bpms[bpm],
            ),
        )
        plot(
            y,
            label="y[b=%d] %s"
            % (
                bunch,
                self.bpms[bpm],
            ),
        )
        title(dumpdate(self.dt))
        legend()

    def plot_xy_fft(self, bpm, bunch):
        x, y = self.get_xy(bpm, bunch)
        f = linspace(0, 0.5, len(x) / 2 + 1)
        plot(
            f,
            abs(rfft(x)),
            label="x[b=%d] %s"
            % (
                bunch,
                self.bpms[bpm],
            ),
        )
        plot(
            f,
            abs(rfft(y)),
            label="y[b=%d] %s"
            % (
                bunch,
                self.bpms[bpm],
            ),
        )
        title(dumpdate(self.dt))
        legend()

    def get_tune(bpm, bunch, tune_fit=maxharm_brent):
        x, y = self.get_xy(bpm, bunch)
        return tune_fit(x - x.mean())[0], tune_fit(y - y.mean())[0]

    def get_tunes(self, tune_fit=maxharm_brent):
        qx = []
        qy = []
        x = self.xbpm
        y = self.ybpm
        for b in range(self.bunches):
            for p in range(self.nbpms):
                xx = x[p, b, :]
                yy = y[p, b, :]
                if sum(xx**2) > 0:
                    qx.append(tune_fit(xx - xx.mean(), 0.20, 0.40))
                if sum(yy**2) > 0:
                    qy.append(tune_fit(yy - yy.mean(), 0.20, 0.40))
        return array(qx), array(qy)

    def get_decay_lin_bunch(self, b):
        taux = []
        tauy = []
        x = self.xbpm
        y = self.ybpm
        for p in range(self.nbpms):
            xx = x[p, b, :]
            yy = y[p, b, :]
            taux.append(get_decay_lin(xx))
            tauy.append(get_decay_lin(yy))
        return array(taux), array(tauy)

    def plot_tune_hist(self):
        tunex, tuney = self.get_tunes()
        qx, ax, px, rx = list(zip(*tunex))
        qy, ay, py, ry = list(zip(*tuney))
        hist(qx, bins=10, label="Qx")
        hist(qy, bins=10, label="Qy")
        legend()
        xlabel("Tune")
        ylabel("Count")
        title(dumpdate(self.dt))
        return self

    def plot_tunes(self, tune_fit=maxharm_brent):
        for b in range(self.bunches):
            tunex, tuney = self.get_tunes_bunch(b)
            qx, ax, px, rx = list(zip(*tunex))
            qy, ay, py, ry = list(zip(*tuney))
            plot(qx, qy, ".")

    def get_coupling(self):
        out = []
        for b in range(self.bunches):
            for p in range(self.nbpms):
                print(b, p)
                xx = self.xbpm[p, b, :]
                yy = self.ybpm[p, b, :]
                cc = fit_coupled_lsq2(xx, yy)
                out.append([b, p] + list(cc))
        return out

    def get_tunes_bunch(self, bunch, tune_fit=maxharm_brent, turns=None):
        qx = []
        qy = []
        x = self.xbpm
        y = self.ybpm
        bpms, bunches, turns = x.shape
        for p in range(self.nbpms):
            xx = x[p, bunch, :]
            yy = y[p, bunch, :]
            if xx[0] == 0 and yy[0] == 0:
                xx = xx[1:]
                yy = yy[1:]
            # if sum(xx**2)>0:
            xx = xx[:turns]
            yy = yy[:turns]
            qx.append(tune_fit(xx - xx.mean()))
            # if sum(yy**2)>0:
            qy.append(tune_fit(yy - yy.mean()))
        return array(qx), array(qy)

    def get_tune_fit(self, tune_fit=maxharm_brent, tunecut=0.03):
        out = []
        for bb in range(self.bunches):
            tunex, tuney = self.get_tunes_bunch(bb, tune_fit=tune_fit)
            out.append(list(zip(*tunex)) + list(zip(*tuney)))
        return out

    def plot_tune_vs_amp(self):
        out = self.get_tune_fit()
        qx = out[:, 0, :].mean(axis=1)
        ax = out[:, 1, :].mean(axis=1)
        qy = out[:, 4, :].mean(axis=1)
        ay = out[:, 5, :].mean(axis=1)
        figure(figsize=(12, 6))
        suptitle("%s" % (self.date))
        subplot(121)
        hexbin(ax, qx, extent=[0.1, 0.4, 0.25, 0.29])
        ylabel("Horizontal tune")
        xlabel("Amplitude [a.u]")
        subplot(122)
        hexbin(ay, qy, extent=[0.1, 0.4, 0.28, 0.32])
        ylabel("Vertical tune")
        xlabel("Amplitude [a.u]")
        return self

    def get_tune_stats(self, tune_fit=maxharm_brent, tunecut=0.03, turns=None):
        qqx = []
        qqy = []
        qqxs = []
        qqys = []
        for bb in range(self.bunches):
            fitx, fity = self.get_tunes_bunch(bb, tune_fit=tune_fit, turns=turns)
            qx, ax, px, rx = fitx.T
            qy, ay, py, ry = fity.T
            if len(qx) > 0:
                qx = qx[abs(qx - 0.28) < tunecut]
                qxm = qx.mean()
                qxs = qx.std()
            else:
                qxm = 0
                qxs = 0
            qqx.append(qxm)
            qqxs.append(qxs)
            if len(qy) > 0:
                qy = qy[abs(qy - 0.31) < tunecut]
                qym = qy.mean()
                qys = qy.std()
            else:
                qxm = 0
                qxs = 0
            qqy.append(qym)
            qqys.append(qys)
        return array([qqx, qqy, qqxs, qqys])

    def get_amp_stats(self, tune_fit=maxharm_brent, tunecut=0.03):
        qqx = []
        qqy = []
        qqxs = []
        qqys = []
        for bb in range(self.bunches):
            fitx, fity = self.get_tunes_bunch(bb, tune_fit=tune_fit)
            qx, ax, px, rx = fitx.T
            qy, ay, py, ry = fity.T
            if len(qx) > 0:
                qx = ax[abs(qx - 0.28) < tunecut]
                qxm = qx.mean()
                qxs = qx.std()
            else:
                qxm = 0
                qxs = 0
            qqx.append(qxm)
            qqxs.append(qxs)
            if len(qy) > 0:
                qy = ay[abs(qy - 0.31) < tunecut]
                qym = qy.mean()
                qys = qy.std()
            else:
                qxm = 0
                qxs = 0
            qqy.append(qym)
            qqys.append(qys)
        return array([qqx, qqy, qqxs, qqys])

    def plot_tune_stats(self, tune_fit=maxharm_brent, tunecut=0.03, turns=None):
        qqx, qqy, qqxs, qqys = self.get_tune_stats(tune_fit, tunecut, turns=turns)
        figure(figsize=(12, 4))
        # suptitle("%s"%(self.date),fontsize=18)
        s1 = subplot(121)
        title("%s" % (self.date))
        index = [j + 1 for j in range(self.bunches)]
        errorbar(index, qqx, yerr=qqxs, label="mean and rms $Q_x$", ecolor="r")
        for xx in [35, 71, 107, 143]:
            axvline(xx, color="g")
        legend(loc=0)
        xlabel("Bunch index")
        ylabel("Horizontal tune")
        y1a, y1b = s1.get_ylim()
        s2 = subplot(122)
        errorbar(index, qqy, yerr=qqys, label="mean and rms $Q_y$", ecolor="r")
        for xx in [35, 71, 107, 143]:
            axvline(xx, color="g")
        legend(loc=0)
        xlabel("Bunch index")
        ylabel("Vertical  tune")
        y2a, y2b = s2.get_ylim()
        s1.set_ylim(min(y1a, y2a), max(y1b, y2b))
        tight_layout()
        return self

    def plot_amp_stats(self, tune_fit=maxharm_brent, tunecut=0.03):
        qqx, qqy, qqxs, qqys = self.get_amp_stats(tune_fit, tunecut)
        figure(figsize=(12, 4))
        # suptitle("%s"%(self.date),fontsize=18)
        s1 = subplot(121)
        title("%s" % (self.date))
        index = [j + 1 for j in range(self.bunches)]
        errorbar(index, qqx, yerr=qqxs, label="mean and rms amplitude", ecolor="r")
        for xx in [35, 71, 107, 143]:
            axvline(xx, color="g")
        legend(loc=0)
        xlabel("Bunch index")
        ylabel("Horizontal amplitude [au]")
        y1a, y1b = s1.get_ylim()
        s2 = subplot(122)
        errorbar(index, qqy, yerr=qqys, label="mean and rms amplitude", ecolor="r")
        for xx in [35, 71, 107, 143]:
            axvline(xx, color="g")
        legend(loc=0)
        xlabel("Bunch index")
        ylabel("Vertical  amplitude [au]")
        y2a, y2b = s2.get_ylim()
        s1.set_ylim(min(y1a, y2a), max(y1b, y2b))
        s2.set_ylim(min(y1a, y2a), max(y1b, y2b))
        draw()
        tight_layout()
        return self
