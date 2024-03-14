from numpy import *
from matplotlib.pyplot import *
from glob import glob
import os

from foot import *

figure(figsize=(4.2, 4))
title("no feed-down")
f2 = FootTrack("b6_100_nocrossing_nocorrection_nosextupoles")
f1 = FootComp("b6_100_nocrossing_nocorrection_nosextupoles")
f2.plot_footprint(nsigma=12, spread=0.015, label="no corr., tracking")
f1.plot_footprint(nsigma=12, spread=0.015, label="no corr., formulas")

f2 = FootTrack("b6_100_nocrossing_correction_nosextupoles")
f1 = FootComp("b6_100_nocrossing_correction_nosextupoles", mcx=True)
f2.plot_footprint(nsigma=12, spread=0.015, label="corr., tracking")
f1.plot_footprint(nsigma=12, spread=0.015, label="corr., formulas")

legend()
tight_layout()
xlim(0.270, 0.290)
ylim(0.300, 0.330)
savefig("footprint_nocross.png")
gca().legend().set_visible(False)
xlim(0.2790, 0.2805)
ylim(0.3092, 0.3109)
xticks(arange(0.2790, 0.2805, 0.0005))
savefig("footprint_nocross2.png")


figure(figsize=(4.2, 4))
title("feed-down")
f2 = FootTrack("b6_100_crossing_nocorrection_nosextupoles")
f1 = FootComp("b6_100_crossing_nocorrection_nosextupoles")
f2.plot_footprint(nsigma=12, spread=0.015, label="no corr., tracking")
f1.plot_footprint(nsigma=12, spread=0.015, label="no corr., formulas")

f2 = FootTrack("b6_100_crossing_correction_nosextupoles")
f1 = FootComp("b6_100_crossing_correction_nosextupoles", mcx=True)
f2.plot_footprint(nsigma=12, spread=0.015, label="corr., tracking")
f1.plot_footprint(nsigma=12, spread=0.015, label="corr., formulas")

# legend()
tight_layout()
xlim(0.270, 0.290)
ylim(0.300, 0.330)
savefig("footprint_cross.png")
# gca().legend().set_visible(False)
xlim(0.2790, 0.2805)
ylim(0.3092, 0.3109)
xticks(arange(0.2790, 0.2805, 0.0005))
savefig("footprint_cross2.png")
