from pymad import *

f = Frame()
f.kq1 = 0.005
f.kq2 = -0.005
f.lq = 4

f.ls = Elem(l=20)
f.q1 = Elem(l=expr("lq"), kn1l=expr("kq1*lq"))
f.q2 = Elem(l=expr("lq*2"), kn1l=expr("kq2*l"))
f.q3 = Elem(l=expr("lq"), kn1l=expr("kq1*lq"))
f.le = Elem(l=4)

f.seq = Line("ls", "q1", "q2", "q3", "le")
f.init = Beam(betx=0.5, bety=0.25)

tt = f.seq.track(f.init)
tt.print_table("name s betx bety alfx alfy")


import scipy.optimize


def ftosolve(x):
    f.kq1, f.kq2 = x
    tt = f.seq.track(f.init)
    return tt.le.alfx, tt.le.alfy


def mk_dir():
    vec = random.rand(2) - 0.5
    return vec / sqrt(sum(vec**2))


x0 = scipy.optimize.fsolve(ftosolve, array([0.005, -0.005]), xtol=1e-13)


out = []
dv = arange(0, 1, 0.01)
for d in dv:
    x1 = x0 + mk_dir() * d
    info = scipy.optimize.fsolve(ftosolve, x1, full_output=1, xtol=1e-13)
    if info[2] == 1:
        nfev = info[1]["nfev"]
        out.append(nfev)
    else:
        out.append(-1)

plot(dv, out)

# 1) find a suitable range for dv to explore the features of the dv,out plot
# 2) check whether options `epsfcn` has any effect
# 3) try to use `scipy.optimize.broyden1` or `scipy.optimize.broyden2` instead of fsolve
