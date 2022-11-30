from numpy import *
import numpy as np


def solveconst(A, b, C, d):
    """Solve constrained least square problem using Lagrange multipliers
    min(|| A x - b ||_2) and Cx=b
    x: n
    A: m x n
    b: m
    C: l x n
    d: l

    L = x A.T A x - 2 A.T x + b.T b + l C x
    Equivalent to:
    (A.T A   C.T )  (x)  = (A.T b)
    (C       0   )  (l)  = (d)
    """
    nl, nx = C.shape
    m = hstack([dot(A.T, A), C.T])
    m = vstack([m, hstack([C, zeros((nl, nl))])])
    n = hstack([dot(A.T, b), d])
    sol = linalg.solve(m, n)
    return sol[:nx]


def makeA(x, N):
    return column_stack([x**i for i in range(N + 1)])


def makeAp(x, N):
    return column_stack([zeros(len(x))] + [i * x ** (i - 1) for i in range(1, N + 1)])


def makeApp(x, N):
    z = [zeros(len(x))] * 2
    z += [i * (i - 1) * x ** (i - 2) for i in range(2, N + 1)]
    return column_stack(z)


def poly_val(p, x):
    return sum([p[i] * x**i for i in range(len(p))], axis=0)


def poly_print(p, x="x", power="**", mul="*"):
    res = ["%.10e" % p[0]]
    if len(p) > 1:
        res.append("%+.10e%s%s" % (p[1], mul, x))
    for i in range(2, len(p)):
        res.append("%+.10e%s%s%s%d" % (p[i], mul, x, power, i))
    return "".join(res)


def poly_fit(N, xdata, ydata, x0=[], y0=[], xp0=[], yp0=[], xpp0=[], ypp0=[]):
    A = makeA(xdata, N)
    b = ydata
    C0 = makeA(array(x0), N)
    C1 = makeAp(array(xp0), N)
    C2 = makeApp(array(xpp0), N)
    C = vstack([C0, C1, C2])
    d = hstack([y0, yp0, ypp0])
    p = solveconst(A, b, C, d)
    return p


if __name__ == "__main__":
    from matplotlib.pyplot import *

    x = linspace(0, 1.0, 101)
    p = array([2, 3, -8, 5])
    y = poly_val(p, x)
    n = 0.02 * sin(2 * pi * 8 * x)
    N = len(p) - 1
    clf()
    plot(x, y + n)
    plot(x, poly_val(poly_fit(4, x, y + n, [0, 1], [2, 2]), x))
    A = makeA(x, N)
    np.linalg.lstsq(A, y + n)
    xx = array([0, 1])
    yy = 2 + 3 * xx - 8 * xx**2 + 5 * xx**3
    C = makeA(xx, N)
    d = yy
    p = solveconst(A, y + n, C, d)
    A = array([[2.0, 3], [1, 2], [4, 5]])
    b = array([2, 3, 1])
    C = array([[1, 0], [0, 1]])
    d = array([2, 3])
    C = array([[1, 0]])
    d = array([2])
    p = solveconst(A, b, C, d)
