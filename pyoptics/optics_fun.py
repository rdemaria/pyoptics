from math import *


def twiss2map(bet1, alf1, bet2, alf2, mu):
    b1b2 = sqrt(bet1 * bet2)
    b1onb2 = sqrt(bet1 / bet2)
    c = cos(2 * pi * mu)
    s = sin(2 * pi * mu)
    r11 = (c + alf1 * s) / b1onb2
    r12 = b1b2 * s
    r21 = ((alf1 - alf2) * c - (1 + alf1 * alf2) * s) / b1b2
    r22 = b1onb2 * (c - alf2 * s)
    return [[r11, r12], [r21, r22]]
