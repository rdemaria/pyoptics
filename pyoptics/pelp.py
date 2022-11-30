from numpy import *


def pelp(time, Ii, If, A, D, R, Te, Ti):
    iaTe = A / 2 * (Te - Ti) ** 2 + Ii
    iapTe = A * (Te - Ti)
    b = iapTe / iaTe
    a = iaTe * exp(-b * Te)
    tl = log(R / (a * b)) / b if Te != 0 else R / A + Ti
    il = a * exp(b * tl) if Te != 0 else A / 2 * (tl - Ti) ** 2 + Ii
    td = (If - il) / R + tl - R / (2 * D)
    tf = td + R / D
    I = (
        A / 2 * (time - Ti) ** 2 + Ii
        if (Ti <= time) and (time < Te)
        else a * exp(b * time)
        if (Te <= time) and (time < tl)
        else R * (time - tl) + il
        if (tl <= time) and (time < td)
        else -D / 2 * (tf - time) ** 2 + If
        if (td <= time) and (time < tf)
        else 0
    )
    return I


if __name__ == "__main__":
    pelp = vectorize(pelp, excluded=[1, 2, 3, 4, 5, 6, 7, 8])
    t = arange(0, 2000)
    p = pelp(t, 750, 11800, 0.009, 0.02, 10.0, 325, 0)
