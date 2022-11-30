from pyoptics import madlang

mad = madlang.open("SPS_Q20_thin.seq")

mad.kmba *= 2
mad.mba_40370__3.knl[0] * 10 == mad.kmba
mad.kmba /= 2
out, rest = mad.sps.expand_struct()
