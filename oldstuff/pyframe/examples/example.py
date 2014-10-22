# Usage:
# python example
# or
# ipython -pylab for interactive use

from pyview import *





# Create general tables

t=table('qd d qf d','x y z')
t.x=1
t.d.y=2
t.qf.all=3
print t


# with nested structure and query routines
tt=table(['start',('f1',t),('f2',t),('f1',t),'end'],'x z')
tt.f1.z=1
tt.show('\.d','x z')

#more complex queries and extraction
tt.show( (tt //'\.d') & (tt.z>0) ,'x z')
print tt[(tt //'\.d') & (tt.z>0)]

#subclass 'tfstable' can load tfs tables:
t1=tfstable('twiss.ir5b1.data.gz')
t1.show(t1//'ip','name ^s$ bet')

print t1.ip5.betx

#the row attributes can be disabled for big tables:
t2=tfstable('twiss.lhcb1.data.gz',rowattr=False)
t2.show(t2//'ip','name ^s$ bet')

print len(t2.s)

#plotting is avaiable via the matplotlib module
from pylab import savefig
t2.plot('betx bety','dx')
savefig('plot.png')

#another table type 'trtable' let compute some basic twiss parameters
# table creation
from trtable import *
s=trtable('d0 q1 d1 q2 d2 q3 d3 q4 d4 end')
# length assigment
s.d0.l, s.d1.l, s.d2.l, s.d3.l, s.d4.l = 23, 2.7, 1, 2.9, 10
s.q1.l, s.q2.l = 9.2, 7.8
# gradient assigment
s.q1.kn1l = 122./23350*s.q1.l
s.q3.l, s.q4.l = s.q2.l, s.q1.l
s.q2.kn1l = -s.q1.kn1l/s.q1.l*s.q2.l
s.q3.kn1l, s.q4.kn1l = s.q2.kn1l, s.q1.kn1l
# initial conditions
s.d0.betx,s.d0.bety=0.25,0.25
# twiss parameter computation
s.track()
# plotting
s.plot('betx bety')
savefig('triplet.png')

#the table can be sliced (overload of * operator)
sl=(s*10)
sl.track()
sl.plot('betx bety')
savefig('triplet.png')

#trtable can be create from a tfs table
s2=from_tfstable(t2,rowattr=False)

# compute  a twiss (30ms on my laptop)
s2.track()

#trtable store the values at the begining of an element
print s2.s[1:]-t2.s
#precision
print '%2f%%' % (abs(1-s2.bety[1:]/t2.bety).max()*100)
print '%2f%%' % (abs(1-s2.betx[1:]/t2.betx).max()*100)




