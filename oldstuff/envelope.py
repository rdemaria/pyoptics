from numpy import zeros,array,sin,cos,dot,sqrt,maximum,average
import pylab as _p

class envelope:
  """
  s1=m.madtable('survey.lhcb1.data')
  s2=m.madtable('survey.lhcb2.data')
  t1=m.madtable('twiss.lhcb1.data')
  t2=m.madtable('twiss.lhcb2.data')
  ap1=m.envelope(s1,t1)
  ap2=m.envelope(s1,t1)

  hold(False)
  figure(figsize=(6,6))
  hold(True)
  plot(ap1.co[:,2],ap1.co[:,0])
  plot(ap1.xp[:,2],ap1.xp[:,0],'g',linewidth=.1)
  plot(ap1.yp2D[:,2],ap1.yp2D[:,0],'r',linewidth=.1)
  plot(ap1.xm[:,2],ap1.xm[:,0],'g',linewidth=.1)
  plot(ap1.ym2D[:,2],ap1.ym2D[:,0],'r',linewidth=.1)
  axis([-10000,10000,-15000,5000])

  t1.select(t1.pattern("IP1"),t1.name,t1.s,ap1.xsize,ap1.ysize)
  savefig('ring.eps',dpi=600)
  """

#    kbeta=1.1,      # beta beating
#    nsigma=9.5,      # 9.5
#    nemit=3.75E-6,   #3.75E-6
#    delta=1.129E-4, # RMS energy spread
#    tol=4.6E-3,  # CO=3mm + dtol=1.6mm
#    deltamax=8E-4,   # for chromaticity measurment
#    betamaxarc=180,  # maximum beta in the arcs
#    dxmaxarc=2,  # maximum beta in the arc
  def __init__(self,s,t, kbeta=1.1, nsigma=10, nemit=3.75E-6, delta=1.129E-4, tol=4.6E-3, deltamax=8.6E-4, betamaxarc=180, dxmaxarc=2, gamma=7400):
    self.co=zeros([len(s.x),3],float)
    self.xp=zeros([len(s.x),3],float)
    self.xm=zeros([len(s.x),3],float)
    self.yp=zeros([len(s.x),3],float)
    self.ym=zeros([len(s.x),3],float)
    self.yp2D=zeros([len(s.x),3],float)
    self.ym2D=zeros([len(s.x),3],float)
    self.xsize=zeros(len(s.x),float)
    self.ysize=zeros(len(s.x),float)

    for i in xrange(len(s.x)):
      vro = array([s.x[i],s.y[i],s.z[i]])
      theta,phi,psi = s.theta[i],s.phi[i],s.psi[i]
      betx,bety,dx,dy,x,y= t.betx[i],t.bety[i],t.dx[i],t.dy[i],t.x[i],t.y[i]
      thetam=array([[cos(theta) ,           0,sin(theta)],
           [          0,           1,         0],
           [-sin(theta),           0,cos(theta)]])
      phim=  array([[          1,          0,          0],
          [          0,cos(phi)   ,   sin(phi)],
          [          0,-sin(phi)  ,   cos(phi)]])
      psim=  array([[   cos(psi),  -sin(psi),          0],
          [   sin(psi),   cos(psi),          0],
          [          0,          0,          1]])
      wm=dot(thetam,dot(phim,psim))
      ex=dot(wm,array([1,0,0]))
      ey=dot(wm,array([0,1,0]))
      self.co[i]=vro+x * ex + y * ey
      emit=nemit/gamma
      dx= dxmaxarc*sqrt(betx/betamaxarc)*0.27+abs(dx)
      dy= dxmaxarc*sqrt(bety/betamaxarc)*0.27+abs(dy)
      xsize=kbeta* (  nsigma*sqrt(betx*emit) + deltamax*dx  ) + tol
      ysize=kbeta* (  nsigma*sqrt(bety*emit) + deltamax*dy  ) + tol
      self.xp[i]=self.co[i] + xsize * ex
      self.xm[i]=self.co[i] - xsize * ex
      self.yp[i]=self.co[i] + ysize * ey
      self.ym[i]=self.co[i] - ysize * ey
      self.yp2D[i]=self.co[i] + ysize * ex
      self.ym2D[i]=self.co[i] - ysize * ex
      self.xsize[i]=xsize
      self.ysize[i]=ysize



def plotenvelopex(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
  plotenvelope("ip5",t1,ap1,"mqy_4l5_b1","mqy_4r5_b1",t2,ap2,"mqy_4l5_b2","mqy_4r5_b2")
  """
  #select
  yip5=(ap1.co[t1._row_ref[ip],0]+ap2.co[t2._row_ref[ip],0])/2
  xip5=(ap1.co[t1._row_ref[ip],2]+ap2.co[t2._row_ref[ip],2])/2
  idxs1=t1._row_ref[sele]
  idxe1=t1._row_ref[eele]
  idxs2=t2._row_ref[sele2]
  idxe2=t2._row_ref[eele2]
  # start plot
  _p.hold(True)
  _p.title("Horizontal beam envelope")
  _p.xlabel(r"$z [\rm{m}]$")
  _p.ylabel(r"$x [\rm{m}]$")
  _p.grid(True)
  # closed orbit
  x1=ap1.co[idxs1:idxe1,2]-xip5
  y1=ap1.co[idxs1:idxe1,0]-yip5
  _p.plot(y1,x1,color=[0,0,1])
  x2=ap2.co[idxs2:idxe2,2]-xip5
  y2=ap2.co[idxs2:idxe2,0]-yip5
  _p.plot(y2,x2,color=[1,0,0])
  # beam1
  x1=ap1.xp[idxs1:idxe1,2]-xip5
  y1=ap1.xp[idxs1:idxe1,0]-yip5
  x2=ap1.xm[idxs1:idxe1,2]-xip5
  y2=ap1.xm[idxs1:idxe1,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='b',alpha=0.2)
  # beam2
  x1=ap2.xp[idxs2:idxe2,2]-xip5
  y1=ap2.xp[idxs2:idxe2,0]-yip5
  x2=ap2.xm[idxs2:idxe2,2]-xip5
  y2=ap2.xm[idxs2:idxe2,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='r',alpha=0.2)


def plotenvelopey(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
  plotenvelopey("ip1",t1,ap1,"mqy_4l1_b1","mqy_4r1_b1",t2,ap2,"mqy_4l1_b2","mqy_4r1_b2")
  """
  #select
  yip5=(ap1.co[t1._row_ref[ip],0]+ap2.co[t2._row_ref[ip],0])/2
  xip5=(ap1.co[t1._row_ref[ip],1]+ap2.co[t2._row_ref[ip],1])/2
  idxs1=t1._row_ref[sele]
  idxe1=t1._row_ref[eele]
  idxs2=t2._row_ref[sele2]
  idxe2=t2._row_ref[eele2]
  # start plot
  _p.hold(True)
  _p.title("Vertical beam envelope")
  _p.xlabel(r"$z [\rm{m}]$")
  _p.ylabel(r"$y [\rm{m}]$")
  _p.grid(True)
  # closed orbit
  x1=ap1.co[idxs1:idxe1,1]-xip5
  y1=ap1.co[idxs1:idxe1,0]-yip5
  _p.plot(y1,x1,color=[0,0,1])
  x2=ap2.co[idxs2:idxe2,1]
  y2=ap2.co[idxs2:idxe2,0]-yip5
  _p.plot(y2,x2,color=[1,0,0])
  # beam1
  x1=ap1.yp[idxs1:idxe1,1]-xip5
  y1=ap1.yp[idxs1:idxe1,0]-yip5
  x2=ap1.ym[idxs1:idxe1,1]-xip5
  y2=ap1.ym[idxs1:idxe1,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='b',alpha=0.2)
  # beam2
  x1=ap2.yp[idxs2:idxe2,1]-xip5
  y1=ap2.yp[idxs2:idxe2,0]-yip5
  x2=ap2.ym[idxs2:idxe2,1]-xip5
  y2=ap2.ym[idxs2:idxe2,0]-yip5
  x = _p.concatenate( (x1,x2[::-1]) )
  y = _p.concatenate( (y1,y2[::-1]) )
  _p.fill(y, x, facecolor='r',alpha=0.2)


def plotbeamsep(s1,t1,s2,t2,ap1,ap2,eps=5E-10):
  idx1=t1// 'bbk_.*[lr]1'
  idx2=t2// 'bbk_.*[lr]1'
# plot(-ap1.co[idx1][:,0],ap1.co[idx1][:,2],'ob')
# plot(-ap2.co[idx2][:,0],ap2.co[idx2][:,2],'or')
  x=ap1.co[idx1][:,0]-ap2.co[idx2][:,0]
  y=ap1.co[idx1][:,1]-ap2.co[idx2][:,1]
  z=ap1.co[idx1][:,2]-ap2.co[idx2][:,2]
  ds1=sqrt(x**2+y**2+z**2)
  bx=maximum(t1.betx[idx1],t2.betx[idx2])
  by=maximum(t1.bety[idx1],t2.bety[idx2])
  s=sqrt(by*eps)
  sep1=ds1/s
  idx1=t1// 'bbk_.*[lr]5'
  idx2=t2// 'bbk_.*[lr]5'
#  plot(-ap1.co[idx1][:,0],ap1.co[idx1][:,2],'ob')
#  plot(-ap2.co[idx2][:,0],ap2.co[idx2][:,2],'or')
  x=ap1.co[idx1][:,0]-ap2.co[idx2][:,0]
  y=ap1.co[idx1][:,1]-ap2.co[idx2][:,1]
  z=ap1.co[idx1][:,2]-ap2.co[idx2][:,2]
  ds5=sqrt(x**2+y**2+z**2)
  bx=maximum(t1.betx[idx1],t2.betx[idx2])
  by=maximum(t1.bety[idx1],t2.bety[idx2])
  s=sqrt(by*eps)
  sep5=ds5/s
  sb=t1.s[idx1]-average(t1.s[idx1])
  _p.plot(sb,sep1,label='ip1')
  _p.plot(sb,sep5,label='ip5')
  _p.title('Beam beam separation')
  _p.xlabel('s [m]')
  _p.ylabel(r'$\sigma$')
  _p.ylim(0,20)
  _p.legend()
  _p.grid()

