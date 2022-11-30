import os
import time
import copy

import matplotlib.pyplot as pl
import matplotlib
import numpy as np


from .optics import *


class Scenarios(object):
    def __init__(self, data):
        self._data = data
        for name in data["scenario_list"]:
            scen = Scenario(name, **data[name])
            setattr(self, name, scen)


class Scenario(object):
    beam_data = ["p", "E", "N<sub>b</sub>", "k<sub>b</sub>", "&epsilon;"]
    beam_data_unit_tmp = ["", "GeV", "10<sup>%d</sup>", "", "&mu;m"]
    ip_data = [
        "&beta;<sub>x</sub>",
        "&beta;<sub>y</sub>",
        "x",
        "y",
        "p<sub>x</sub>",
        "p<sub>y</sub>",
    ]
    ip_data_unit = ["m", "m", "mm", "mm", "&mu;rad", "&mu;rad"]
    ip_names = ["ip1", "ip2", "ip5", "ip8"]

    def __init__(self, name, **data):
        self.name = name
        self.__dict__.update(data)
        self.beam_data_unit = Scenario.beam_data_unit_tmp[:]
        self.beam_data_unit[2] = self.beam_data_unit_tmp[2] % (self.npart_unit)
        for cname, cdata in list(self.confs.items()):
            Configuration._instances[(name, cname)] = cdata
        for cname, cdata in list(self.confs.items()):
            self.confs[cname] = Configuration(
                name=cname, scenario=name, scn=self, **cdata
            )
            setattr(self, cname, self.confs[cname])


class Configuration(object):
    selection = [
        "LHC",
        "IR1",
        "Arc12",
        "IR2",
        "Arc23",
        "IR3",
        "Arc34",
        "IR4",
        "Arc45",
        "IR5",
        "Arc56",
        "IR6",
        "Arc67",
        "IR7",
        "Arc78",
        "IR8",
        "Arc81",
    ]
    _instances = {}
    pdata = {"p": (0.931494061, 1), "Pb": (193.68715, 82)}  # GeV,charge  # GeV,charge

    def __init__(self, **data):
        scenario = data["scenario"]
        if "template" in data:
            tname = (scenario, data["template"])
            tmp = self._instances[tname].copy()
            self.__dict__.update(copy.deepcopy(tmp))
        self.__dict__.update(data)
        if hasattr(self, "settings"):
            part1, self.nrj1, self.np1, self.nb1, self.emit_n1 = self.settings["beam1"]
            part2, self.nrj2, self.np2, self.nb2, self.emit_n2 = self.settings["beam2"]
            if part1 == "p":
                self.part1 = "proton"
            else:
                self.part1 = "ion"
            if part2 == "p":
                self.part2 = "proton"
            else:
                self.part2 = "ion"
            self.part1_web = part1
            self.part2_web = part2
            self.pmass1, self.charge1 = self.pdata[part1]
            self.pmass2, self.charge2 = self.pdata[part2]
            self.emit1 = self.emit_n1 / self.nrj1 * self.pmass1
            self.emit2 = self.emit_n2 / self.nrj2 * self.pmass2
            self.np1_web = float(self.np1) / 10**self.scn.npart_unit
            self.np2_web = float(self.np2) / 10**self.scn.npart_unit
            self.emit_n1_web = self.emit_n1 * 1e6
            self.emit_n2_web = self.emit_n2 * 1e6

    def get_conf_dir(self):
        return os.path.join(basedir, self.scenario, self.name)

    def get_twiss_fn(self, beam):
        return os.path.join(self.get_conf_dir(), "twiss_lhcb%s.tfs" % beam)

    def get_twiss(self, beam):
        return optics.open(self.get_twiss_fn(beam))

    def get_data(self):
        for b12 in "12":
            t = optics.open(self.get_twiss_fn(b12))
            for nn, ipn in enumerate([1, 2, 5, 8]):
                print((get_ip_data(t, ipn)))
                # exp=self.settings['exp'][nn]
                # self.settings['ip%sb%s'%(ipn,b12)]=get_ip_data(t,ipn)+[exp]
                self.settings["ip%sb%s" % (ipn, b12)] = get_ip_data(t, ipn)
        return self


def get_ip_data(t, n):
    i = np.where(t // ("IP%d" % n))[0][0]
    data = (
        t.betx[i],
        t.bety[i],
        t.x[i] * 1e3,
        t.y[i] * 1e3,
        t.px[i] * 1e6,
        t.py[i] * 1e6,
    )
    # print data
    return [s2d_conv(i) for i in data]


def s2d_conv(n):
    if abs(n) > 1e-6:
        n = round(n, -int(np.log10(abs(n))) + 2)
    else:
        n = 0
    s = ("%.3f" % n).replace(".000", "")
    if s == "-0":
        s = "0"
    return s


LHCopt = None
basedir = "/afs/cern.ch/work/l/lhcopt/public/lhc_optics_web/www"
datafn = os.path.join(basedir, "data.json")
if os.path.exists(datafn):
    # try:
    import json

    data = json.load(open(datafn))
    LHCopt = Scenarios(data)
# except:
#  pass
