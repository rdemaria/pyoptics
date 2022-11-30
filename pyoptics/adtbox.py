import time

from numpy import *
import scipy.signal
from matplotlib.pyplot import *
import h5py

from pyoptics import harmonic_fit_2


def fft2d(vv, fv):
    ff = empty((len(fv), len(vv.T)))
    for iii, vvv in enumerate(vv.T):
        zz = where(diff(vvv == 0) != 0)[0]
        if len(zz) > 0:
            vvv = vvv[zz[0] :]
            # print zz[0]
        ffv = linspace(0, 0.5, len(vvv) / 2 + 1)
        ff[:, iii] = interp(fv, ffv, abs(fft.rfft(vvv * 1.0)))
    return ff


class ADTBox(object):
    def __init__(self, filename):
        fh = h5py.File(filename, "r")
        self.filename = filename
        if "B1" in fh:
            self.beam = 1
        else:
            self.beam = 2
        self.ah = fh["B%d/horizontal" % self.beam][:]
        self.av = fh["B%d/vertical" % self.beam][:]
        self.label = filename.split("/")[-1].replace("ADTObsBox_", "")[:-3]
        tm = time.strptime(self.label, "%Y%m%d_%H%M%S_%f")
        self.timestamp = time.mktime(tm) + float(self.label.split("_")[2]) / 1e6

    def get_fft2d(self, hv, fv):
        if hv == "h":
            vv = self.ah
        elif hv == "v":
            vv = self.av
        ff = fft2d(vv, fv)
        return ff

    def plot_spectra_2d(self):
        sub = 121
        suptitle("Beam %d: %s" % (self.beam, self.label))
        fv = linspace(0, 0.5, self.ah.shape[0] / 2 + 1)
        out = []
        for hv in "hv":
            subplot(sub)
            ff = self.get_fft2d(hv, fv)
            imshow(abs(ff), aspect="auto", origin="bottom", extent=[0, 3564, 0, 0.5])
            mm = abs(ff[abs(fv - 0.3) < 0.1]).max()
            clim(0, mm)
            sub += 1
            xlabel("slot")
            ylabel(r"freq %s [$\rm f_{rev}$]" % hv.upper())
            out.append(ff)
        return out

    def plot_osc_2d(self):
        sub = 121
        suptitle("Beam %d: %s" % (self.beam, self.label))
        fv = linspace(0, 0.5, self.ah.shape[0] / 2 + 1)
        for hv in "hv":
            vv = getattr(self, "a" + hv)
            subplot(sub)
            imshow(vv, aspect="auto", origin="bottom")
            sub += 1
            xlabel("slot")
            ylabel("turns %s" % hv.upper())

    def plot_env_2d(self, na, nb):
        sub = 121
        suptitle("Beam %d: %s" % (self.beam, self.label))
        fv = linspace(0, 0.5, self.ah.shape[0] / 2 + 1)
        out = []
        for hv in "hv":
            vv = getattr(self, "a" + hv)[na:nb].copy()
            vv -= vv.mean(axis=0)
            hb = abs(scipy.signal.hilbert(vv, axis=0))
            subplot(sub)
            imshow(hb, aspect="auto", origin="bottom")
            sub += 1
            xlabel("slot")
            ylabel("turns %s" % hv.upper())
            out.append(hb)
        return out
