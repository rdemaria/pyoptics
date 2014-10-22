
em,nb,c3,inj,c5,rmp,sq=loadtxt('IBS.dat').reshape(100,100,7).transpose(2,0,1)

em*=1e6;
nb/=1e11;
eg=(1+inj)*(1+rmp)*(1+sq)




#ext=[em[0,0],em[-1,0],nb[0,0],nb[0,-1]]
#imshow(eg*em,origin='top',aspect='auto',extent=ext)

#ylabel(r'$\sqrt{\epsilon_x \epsilon_y}/\gamma {\rm [\mu m]}$ ')

cs=contourf(nb,em,eg*em,arange(1,2.8,0.2))
ylabel(r'$\epsilon_{\rm inj}/\gamma {\rm [\mu m]}$ ')
xlabel(r'$N_b [10^{11}]$ ')
cb=colorbar()
cb.set_label(r'$\epsilon_{\rm col}/\gamma {\rm [\mu m]}$')

title(r'$\epsilon$ blow-up after injection, ramp and squeeze ($\sigma_t=10$ cm)')


def fm(x):
  br=nb/em
  #return em+x[0]*nb*x[1]/em**x[2]+
  return em+x[0]*nb**x[1]/em**x[2]

def ftomin(x):
  return sum((eg*em-fm(x))**2)

import scipy.optimize

x0=[0.2,1,1]
x1=scipy.optimize.fmin(ftomin,x0)
print x1

clf()
contour(nb,em,eg*em,arange(1,2.8,0.2))
contour(nb,em,fm(x1),arange(1,2.8,0.2))

title(r'$\epsilon_{\rm col}=\epsilon_{\rm inj}+%.1f \, N_b^{%.1f} / \epsilon_{\rm inj}^{%.1f}$'%tuple(x1))


