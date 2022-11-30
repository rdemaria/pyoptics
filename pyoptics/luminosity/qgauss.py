def gauss(x, sig):
    return 1 / sqrt(2 * pi) / sig * exp(-(x**2) / (2 * sig**2))


def qgauss(x, sig):
    S = sig * 4 * sqrt(2)
    return 32 / (5 * pi * S) * (1 - 4 * (x / S) ** 2) ** 2.5
