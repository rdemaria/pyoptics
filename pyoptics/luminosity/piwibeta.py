from numpy import *
import scipy
from matplotlib.pyplot import *

from pyoptics.luminosity.archive.lumi import *


def piwibeta(sigma_z=0.0755, bety=0.075, dsep=10):
    def ftomin(x):
        return -luminosity(
            betx=x[0], bety=bety, sigma_z=sigma_z, dsep=dsep, debug=False
        )

    (betx,) = scipy.optimize.fmin(ftomin, [3])
    return betx


def plp(sigma_z=0.0755, dsep=10):
    x = arange(0.05, 1.0, 0.01)
    y = [piwibeta(bety=i, sigma_z=sigma_z, dsep=dsep) for i in x]
    lbl = r"$\sigma_z=%g$ m, $d_{\rm sep}=%g \sigma$" % (sigma_z, dsep)
    plot(x, y, label=lbl)


plp(sigma_z=0.05, dsep=10)
plp(sigma_z=0.075, dsep=10)
plp(sigma_z=0.10, dsep=10)

plp(sigma_z=0.05, dsep=12)
plp(sigma_z=0.075, dsep=12)
plp(sigma_z=0.10, dsep=12)

plp(sigma_z=0.05, dsep=14)
plp(sigma_z=0.075, dsep=14)
plp(sigma_z=0.10, dsep=14)


xlabel(r"$\beta_y [m]$")
ylabel(r"$\beta_x [m]$")
grid(True)
legend(loc="lower right")
savefig("piwinski_beta.pdf")
savefig("piwinski_beta.png")

# plp(sigma_z=0.0755,alpha=6.5)
# plp(sigma_z=0.07  ,alpha=6.5)
# plp(sigma_z=0.06  ,alpha=6.5)
# plp(sigma_z=0.05  ,alpha=6.5)
#
#
# xlabel(r'$\beta_y [m]$')
# ylabel(r'$\beta_x [m]$')
# grid(True)
# legend(loc='lower right')
#
# savefig('piwinski_beta_13.pdf')
