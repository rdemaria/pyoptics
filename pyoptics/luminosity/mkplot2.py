from numpy import *
import matplotlib.pyplot as plt

mpl.style.use("classic")

from pyoptics.luminosity import pileup


b = pileup.Beam()
plt.figure(figsize=(5, 4))
plt.clf()
b.plot_lumi_int_virt(17, 2.2e11)
plt.savefig("lumi_int_virt17_2.2.png")
plt.clf()
b.plot_lumi_int_virt(12, 2.2e11, 1960, 3.6)
plt.savefig("lumi_int_virt12_2.2.png")
