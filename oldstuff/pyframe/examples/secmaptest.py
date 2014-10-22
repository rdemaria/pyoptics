import numpy
from qdmad import *
import os


numpy.set_printoptions(precision=5)

os.system('madx <sectormap.madx')
t=tfstable('twiss.data')

q=trtable('mb mb dr end')
q.mb.l, q.mb.knl0, q.mb.knl1, q.dr.l =1.2, .9, -.0, 100-2.4
q.track()

t.show(cols='^s$ ^l$')
q.show(cols='^s$ ^l$')

t.show(cols='betx alfx bety  alfy')
q.show(cols='betx alfx bety  alfy')
sum((t.bety[:-1]-q.bety)**2)


t.show(cols='mux muy')

t.show(cols='dx dpx dy dpy')
q.show(cols='^s$ ^l$')

os.system('madx <sectormap.madx')
numpy.set_printoptions(precision=5)
sm=sectormap('sec.data')
print sm.r['mba1']
print transport.quad([1.2,.9,0,0,0,1,7459.])

