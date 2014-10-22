from pyviews import *
from numpy import *

tf=tfstable('twiss.ir4b1.data')
tm=from_tfstable(tf).track()
tr=from_tfstable(tf).track()


quads="mqt_13l4_b1 mqt_12l4_b1 mqtli_11l4_b1 mqml_10l4_b1 mqm_9l4_b1 mqml_8l4_b1 mqm_7l4_b1 mqy_6l4_b1 mqy_5l4_b1 mqy_5r4_b1 mqy_6r4_b1 mqm_7r4_b1 mqml_8r4_b1 mqm_9r4_b1 mqml_10r4_b1 mqtli_11r4_b1 mq_12r4_b1 mq_13r4_b1"

pre= "mqmc_9l4_b1.kn1l=mqm_9l4_b1.kn1l/mqm_9l4_b1.l*mqmc_9l4_b1.l;"
pre+="mqmc_9r4_b1.kn1l=mqm_9r4_b1.kn1l/mqm_9r4_b1.l*mqmc_9r4_b1.l"
post= "ip4, dx dpx; end, betx bety alfx alfy dx dpx mux muy"

x=tr._all[tr.getaddr(quads+',kn1l')]/tr._all[tr.getaddr(quads+',l')]

out=['def mkfs(x):']
for i,n in enumerate(quads.split()):
  out.append('%-21s=x[%d]*tr.%s.l' %('tr.%s.kn1l'%n,i,n))

out.append('tr.mqmc_9l4_b1.kn1l= tr.mqm_9l4_b1.kn1l/tr.mqm_9l4_b1.l*tr.mqmc_9l4_b1.l')
out.append('tr.mqmc_9r4_b1.kn1l= tr.mqm_9r4_b1.kn1l/tr.mqm_9r4_b1.l*tr.mqmc_9r4_b1.l')
idp=tr.getaddr(post)
out.append('tr.track()')
out.append('return array([')
for r,c,v in tr.showaddr(idp):
#  out.append('  tr.%s.%s%+12.5e,' % (r,c,-v))
  out.append('  tr.%s.%s-tm.%s.%s,' % (r,c,r,c))

for i,n in enumerate(quads.split()):
  out.append('  %12.5e-x[%d],' % (220./23350,i))
## alternating focusing
#for i in quads.split():
#  bx,by='tr.%s.betx'%i,'tr.%s.bety'%i
#  if tr[i].betx>tr[i].bety:
#    out.append('  %-21s-%s,' % (bx,by))
#  else:
#    out.append('  %-21s-%s,' % (by,bx))
# alternating focusing
for i in quads.split():
  bx,by='tr.%s.betx'%i,'tr.%s.bety'%i
  out.append('  500-%-21s,' % (bx))
  out.append('  500-%-21s,' % (by))

out.append('])')

#print '\n  '.join(out)

exec '\n  '.join(out)
mkfs.__doc__='\n  '.join(out)
y=mkfs(x)
mask=-asarray(y,dtype=bool)
