from numpy import *
import matplotlib.pyplot as pl

import pileup


b = pileup.Beam()
pl.figure(figsize=(5, 4))
pl.clf()
b.plot_lumi_int_lev()
pl.savefig("lumi_int_lev.png")
pl.clf()
b.plot_lumi_int_virt(20, 2.2e11)
pl.savefig("lumi_int_virt20_2.2.png")
pl.clf()
b.plot_lumi_int_virt(20, 1.4e11)
pl.savefig("lumi_int_virt20_1.2.png")
pl.clf()
b.plot_lumi_int_virt(6, 1.9e11)
pl.savefig("lumi_int_virt6_1.9.png")
pl.clf()
b.plot_lumi_int_virt(18, 1.9e11)
pl.savefig("lumi_int_virt18_1.9.png")
pl.clf()
b.plot_lumi_int_virt(22, 1.9e11)
pl.savefig("lumi_int_virt22_1.9.png")
pl.clf()
b.plot_lumi_int_virt(4, 1.4e11)
pl.savefig("lumi_int_virt4_1.4.png")


b = pileup.Beam()
ppb = linspace(0.0, 3, 100)

pl.figure(figsize=(6, 5))
b.plot_lumi_int_emit(5.06, 3)
pl.plot(ppb, ppb * 0.95)
pl.ylim(0.7, 3)
pl.plot([2, 2], [0, 3])
pl.tight_layout()
pl.savefig("hlstd_emit_ppb.png")

pl.figure(figsize=(6, 5))
b.n_bunches = 2592
b.plot_lumi_int_emit(4.8, 3, cb=1)
pl.plot(ppb, ppb * 0.81 - 0.3)
pl.plot([2, 2], [0, 3])
pl.ylim(0.7, 3)
pl.savefig("hlbcms_emit_ppb.png")
pl.tight_layout()


pl.figure(figsize=(6, 5))
b.plot_lumi_int_emit(1.78, 3, base=20 / 3.6)
pl.plot(ppb, ppb * 1.9)
pl.ylim(0.7, 3)
pl.plot([1.3, 1.3], [0, 3])
pl.tight_layout()
pl.savefig("lhstd_emit_ppb.png")

pl.figure(figsize=(6, 5))
b.n_bunches = 2592
b.plot_lumi_int_emit(1.69, 3, base=20 / 3.6)
pl.plot(ppb, ppb * 1.33 - 0.463)
pl.plot([1.3, 1.3], [0, 3])
pl.ylim(0.7, 3)
pl.tight_layout()
pl.savefig("lhbcms_emit_ppb.png")
