from pyoptics import *

apn='q4_aperture'
apfiles={apn.upper():zip(*[ map(float,l.split()) for l in open(apn)])}

a=BeamEnvelope(optics.open('ap1.tfs'),None,None,apfiles=apfiles)

n1=a.ap.n1[1]
res=a.get_halo_min(1,astep=10)
n2,t1,x0,y0,xp,yp,df,ang,sig=res
print n1,n2,n1-n2


a=BeamEnvelope(optics.open('ap2.tfs'),None,None,apfiles=apfiles)

n1=a.ap.n1[1]
res=a.get_halo_min(1,astep=5)
n2,t1,x0,y0,xp,yp,df,ang,sig=res
print n1,n2,(n1-n2)/n2,ang

a.get_min_dist(1,n2)
a.get_halo_min_m(1,n2)

n3=scipy.optimize.fsolve(lambda x: a.get_min_dist(1,x),5)




clf()
a.plot_aperture(1)
a.plot_pos_tol(1,lbl='tol')
a.plot_halo(1,n1,n1,n1,'r','mad')
a.plot_halo(1,n2,n2,n2,'b','py')
legend(loc=0)
