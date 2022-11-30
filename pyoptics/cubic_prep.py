def cubic_fit1(v0, v1, v0p, v1p, t1):
    a0 = v0
    a1 = v0p
    t2 = t1 * t1
    t3 = t2 * t1
    a2 = (-t1 * (2 * v0p + v1p) - 3 * v0 + 3 * v1) / t2
    a3 = (t1 * (v0p + v1p) + 2 * v0 - 2 * v1) / t3
    return a0, a1, a2, a3


def cubic_fit2(v0, v1, v2, v0p, t1, t2):
    a0 = v0
    a1 = v0p
    det = t1**2 * t2**2 * (t1 - t2)
    det2 = t1**2 - t2**2
    a2 = (
        t1**3 * v2 - t1 * t2 * v0p * det2 - t2**3 * v1 - v0 * (t1**3 - t2**3)
    ) / det
    a3 = (-(t1**2) * v2 + t1 * t2 * v0p * (t1 - t2) + t2**2 * v1 + v0 * det2) / det
    return a0, a1, a2, a3


def cubic_val(t, a0, a1, a2, a3):
    v0 = a0 + t * (a1 + t * (a2 + t * a3))
    v1 = a1 + t * (2 * a2 + 3 * t * a3)
    v2 = 2 * a2 + 6 * t * a3
    return v0, v1, v2


def test_curve():
    from numpy import linspace
    from matplotlib.pyplot import plot

    i = [2.0, 3.0, 4.0, 5.0]
    tt = [0, 2.0, 3.0, 4.0]
    c1 = cubic_fit2(i[0], i[1], i[2], 0, tt[1], tt[2])
    t = linspace(0, 3.0, 100)
    v0, v1, v2 = cubic_val(t, *c1)
    plot(t, v0)
    c2 = cubic_fit1(i[2], i[3], v1[-1], 0, tt[3] - tt[2])
    t = linspace(0, 1.0, 100)
    v0, v1, v2 = cubic_val(t, *c2)
    plot(t + tt[2], v0)
    plot(tt, i, "o")


def get_formulas():
    from sympy import var, solve

    (t, a0, a1, a2, a3) = var("t,a0,a1,a2,a3", real=True)
    (t1, t2, v0, v1, v2, v0p, v1p) = var("t1,t2,v0,v1,v2,v0p,v1p", real=True)

    v = a0 + a1 * t + a2 * t**2 + a3 * t**3
    vp = p.diff(t)
    vv = [a0, a1, a2, a3]
    sol = solve(
        [
            v.subs(t, 0) - v0,
            v.subs(t, t1) - v1,
            vp.subs(t, 0) - v0p,
            vp.subs(t, t1) - v1p,
        ],
        vv,
    )
    for vn in vv:
        print("%s=%s" % (vn, sol[vn]))
    sol = solve(
        [
            v.subs(t, 0) - v0,
            v.subs(t, t1) - v1,
            v.subs(t, t2) - v2,
            vp.subs(t, 0) - v0p,
        ],
        vv,
    )
    for vn in vv:
        print("%s=%s" % (vn, sol[vn]))
