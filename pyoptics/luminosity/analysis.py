from beam import *


b=Beam(N=2e11,alpha=12.27/2,nb=2808,betx=0.15,bety=0.15,emit_n=2.5e-6,sigma_z=0.075)

b.luminosity()


b=Beam(N=2.5e11,alpha=6.5,nb=2808,betx=0.55,bety=0.55,emit_n=3.75e-6,sigma_z=0.0755)
b.luminosity(debug=True)
b.lumi_solve(5e38,'N')
print b.N

