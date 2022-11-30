from __future__ import division

import numpy as np


def quad(l, kn0l, ks0l, kn1l, ks1l, br=1, gr=7000):
    r = np.zeros((6, 6))
    fx = 0
    if l == 0:
        r[0, 0] = 1
        r[0, 1] = 0
        r[1, 0] = -kn1l - kn0l**2
        r[1, 1] = 1
        r[2, 2] = 1
        r[2, 3] = 0
        r[3, 2] = kn1l - ks0l**2
        r[3, 3] = 1
        r[0, 5] = 0
        r[1, 5] = kn0l
    else:
        kx = +kn1l / l + (kn0l / l) ** 2
        if abs(kx) < 1e-10:
            cx = 1 - l**2 * kx / 2
            sx = l - l**3 * kx / 6
            dx = l**2 / 2
            fx = l**3 / 6
        elif kx > 0:
            skx = np.sqrt(kx)
            cx = np.cos(skx * l)
            sx = np.sin(skx * l) / skx
            dx = (1 - cx) / kx
            fx = (l - sx) / kx
        else:
            skx = np.sqrt(-kx)
            cx = np.cosh(skx * l)
            sx = np.sinh(skx * l) / skx
            dx = (1 - cx) / kx
            fx = (l - sx) / kx
        ky = -kn1l / l + (ks0l / l) ** 2
        if abs(ky) < 1e-10:
            cy = 1 - l**2 * ky / 2
            sy = l - l**3 * ky / 6
            dy = l**2 / 2
        elif ky > 0:
            sky = np.sqrt(ky)
            cy = np.cos(sky * l)
            sy = np.sin(sky * l) / sky
            dy = (1 - cy) / ky
        else:
            sky = np.sqrt(-ky)
            cy = np.cosh(sky * l)
            sy = np.sinh(sky * l) / sky
            dy = (1 - cy) / ky
        r[0, 0] = cx
        r[0, 1] = sx
        r[1, 0] = -kx * sx
        r[1, 1] = cx
        r[2, 2] = cy
        r[2, 3] = sy
        r[3, 2] = -ky * sy
        r[3, 3] = cy
        r[0, 5] = dx * kn0l / l
        r[1, 5] = sx * kn0l / l
        r[4, 0] = -r[1, 5]
        r[4, 1] = -r[0, 5]
        r[2, 5] = dy * ks0l / l
        r[3, 5] = sy * ks0l / l
        r[4, 2] = -r[3, 5]
        r[4, 3] = -r[2, 5]
        r[4, 4] = 1
        r[4, 5] = l / br**1 / gr**2 - (kn0l / l) ** 2 * fx / br**2
        r[5, 5] = 1
    return r


def propbeta(r, bet):
    betx, alfx, mux, bety, alfy, muy = bet
    t1 = r[0, 0] * betx - r[0, 1] * alfx
    t2 = r[1, 0] * betx - r[1, 1] * alfx
    newbetx = (t1**2 + r[0, 1] ** 2) / betx
    newalfx = -(t1 * t2 + r[0, 1] * r[1, 1]) / betx
    newmux = mux + np.arctan2(r[0, 1], t1) / (2 * np.pi)
    t3 = r[2, 2] * bety - r[2, 3] * alfy
    t4 = r[3, 2] * bety - r[3, 3] * alfy
    newbety = (t3**2 + r[2, 3] ** 2) / bety
    newalfy = -(t3 * t4 + r[2, 3] * r[3, 3]) / bety
    newmuy = muy + np.arctan2(r[2, 3], t2) / (2 * np.pi)
    return [newbetx, newalfx, newmux, newbety, newalfy, newmuy]


def track(seq, init, br=1, gr=1000):
    out = [["start"] + init]
    # disp=np.zeros((6,1))
    for name, l, kn0l, ks0l, kn1l, ks1l in seq:
        r = quad(l, kn0l, ks0l, kn1l, ks1l, br, gr)
        init = out[-1]
        s = init[1] + l
        bet = init[2:8]
        res = propbeta(r, bet)
        # disp.T[:]=[[dx,dpx,dy,dpy,ct,pt]]
        disp = np.array([init[8:14]]).T
        newdisp = np.dot(r, disp)
        out.append([name, s] + res + list(newdisp.T[0]))
    return out
