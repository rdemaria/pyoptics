# Usage:
# python example
# or
# ipython -pylab for interactive use

from __init__ import *

t2=tfstable('twiss.lhcb1.data.gz',rowattr=False)

s2=from_tfstable(t2,rowattr=False)
s2.track()

#import profile,pstats

#profile.run('[s2.track() for i in xrange(5)]')

import cProfile
cProfile.run('[s2.track() for i in xrange(10)]')

