from pyoptics import *
import aperture
import matplotlib.pyplot as plt


def plot_ap_all(optpath="temp", optmodel="nominal", refap=7):
    aperture.plot_ap(optpath + "/" + optmodel + "/ap_ir1b1.tfs", ref=refap)
    #   plt.savefig(optpath+'/'+optmodel+'/ap1b1.png')
    aperture.plot_ap(optpath + "/" + optmodel + "/ap_ir1b2.tfs", ref=refap)
    #   plt.savefig(optpath+'/'+optmodel+'/ap1b2.png')
    aperture.plot_ap(optpath + "/" + optmodel + "/ap_ir5b1.tfs", ref=refap)
    #   plt.savefig(optpath+'/'+optmodel+'/ap5b1.png')
    aperture.plot_ap(optpath + "/" + optmodel + "/ap_ir5b2.tfs", ref=refap)


#   plt.savefig(optpath+'/'+optmodel+'/ap5b2.png')
