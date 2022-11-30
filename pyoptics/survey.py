import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from . import tfsdata


"""
Primitives: Point, Line, Arc, Text, Curve
Canvas: plot primitives in different projection
Objects: return primimitives
Frame: contain a reference and a list of objects

todo:
- canvas2D
- rethink about render and callbacks

- plot arcs
-  _point to return a Point
- use style dictionary and style names

"""


def _point(point):
    if isinstance(point, Point):
        return point.p
    else:
        return point


class Point:
    def __init__(self, p, color=None, label=None, marker="o"):
        self.p = _point(p)
        self.color = color
        self.label = label
        self.marker = marker

    def render(self):
        return self


class Line:
    def __init__(
        self, p1, p2, label=None, color=None, tips=None, linestyle=None, marker=None
    ):
        self.p1 = _point(p1)
        self.p2 = _point(p2)
        self.label = label
        self.color = color
        self.tips = tips
        self.linestyle = linestyle
        self.marker = marker

    def render(self):
        return self


class Curve:
    def __init__(
        self, xyz, label=None, color=None, tips=None, linestyle=None, marker=None
    ):
        """
        xyz: 3xn array
        """
        self.xyz = xyz
        self.label = label
        self.color = color
        self.tips = tips
        self.linestyle = linestyle
        self.marker = marker

    def render(self):
        return self


class Arc:
    @classmethod
    def from_start_end_center(cls, start, end, center):
        start = _point(start)
        end = _point(end)
        center = _point(center)
        startv = start - center
        endv = end - center
        r1 = normv(startv)
        r2 = normv(endv)
        xv = startv / r1
        ev = endv / r2
        angle = dotp(xv, ev)
        nv = crossp(xv, ev)
        yv = crossp(nv, xv)
        return cls(center, r1 * xv, r2 * yv, angle)

    def __init__(self, center, xv, yv, angle, label=None, linestyle=None, tips=None):
        self.center = _point(center)
        self.xv = _point(xv)
        self.yv = _point(yv)
        self.angle = angle
        self.color = color
        self.label = label
        self.tips = tips
        self.linestyle = linestyle

    def render(self):
        return self


class Text:
    def __init__(self, text, center, color=None):
        self.text = text
        self.center = _point(center)
        self.color = color

    def render(self):
        return self


class Canvas3D:
    """
    labels: labels of the 3D points
    scales: scales of the 3D points
    xyz: mapping index point to index drawing
    """

    @classmethod
    def mad(cls):
        return cls(labels=["X", "Y", "Z"], scale=[-1, 1, 1], xyz=[0, 2, 1])

    def __init__(self, ax=None, labels=["X", "Y", "Z"], scale=[1, 1, 1], xyz=[0, 1, 2]):
        self.labels = labels
        self.xyz = xyz
        self.scale = scale
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            self.ax = ax
            ax.set_xlabel(self.labels[self.xyz[0]])
            ax.set_ylabel(self.labels[self.xyz[1]])
            ax.set_zlabel(self.labels[self.xyz[2]])
            self.set_lim(1)
        self.ax = ax
        self.objects = set()

    def _set_lim(self, a, b, index):
        if self.scale[index] < 0:
            x1, x2 = b, a
        else:
            x1, x2 = a, b
        ff = self.ax.set_xlim, self.ax.set_ylim, self.ax.set_zlim
        ff[self.xyz[index]](x1, x2)

    def set_xlim(self, x1=None, x2=None):
        self._set_lim(x1, x2, 0)

    def set_ylim(self, y1=None, y2=None):
        self._set_lim(y1, y2, 1)

    def set_zlim(self, z1=None, z2=None):
        self._set_lim(z1, z2, 2)

    def set_lim(self, lim):
        self.set_xlim(-lim, lim)
        self.set_ylim(-lim, lim)
        self.set_zlim(-lim, lim)

    def get_xyz(self, p):
        x = p[self.xyz[0]] * np.abs(self.scale[0])
        y = p[self.xyz[1]] * np.abs(self.scale[1])
        z = p[self.xyz[2]] * np.abs(self.scale[2])
        return x, y, z

    def _draw(self, *primitives):
        for pp in primitives:
            childs = pp.render()
            if childs is pp:
                getattr(self, "draw_" + pp.__class__.__name__)(pp)
            else:
                self.draw(*childs)

    def draw(self, *primitives):
        for pp in primitives:
            self.objects.add(pp)
        self._draw(*primitives)
        plt.show()

    def redraw(self):
        self.ax.clear()
        self._draw(*self.objects)

    def draw_Point(self, point):
        x, y, z = self.get_xyz(point.p)
        self.ax.plot([x], [y], [z], color=point.color, marker=point.marker)

    def draw_Line(self, line):
        x1, y1, z1 = self.get_xyz(line.p1)
        x2, y2, z2 = self.get_xyz(line.p2)
        self.ax.plot(
            [x1, x2],
            [y1, y2],
            [z1, z2],
            color=line.color,
            marker=line.marker,
            linestyle=line.linestyle,
        )

    def draw_Text(self, text):
        x, y, z = self.get_xyz(text.center)
        self.ax.text(x, y, z, text.text)

    def draw_Curve(self, curve):
        x, y, z = self.get_xyz(curve.xyz)
        self.ax.plot(
            x, y, z, color=curve.color, marker=curve.marker, linestyle=curve.linestyle
        )


class Canvas2D:
    """
    labels: labels of the 3D points
    scales: scales of the 3D points
    xyz: mapping index point to index drawing
    """

    @classmethod
    def madzx(cls, labels=["X", "Y", "Z"], scale=[1, 1, 1]):
        return cls(labels=labels, scale=scale, xy=[2, 0])

    def __init__(self, ax=None, labels=["X", "Y"], scale=[1, 1, 1], xy=[0, 1]):
        self.labels = labels
        self.xy = xy
        self.scale = scale
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            self.ax = ax
            ax.set_xlabel(self.labels[self.xy[0]])
            ax.set_ylabel(self.labels[self.xy[1]])
            self.set_lim(1)
        self.ax = ax
        self.objects = set()

    def _set_lim(self, a, b, index):
        scale = self.scale[index]
        if index in self.xy:
            if scale < 0:
                x1, x2 = b, a
            else:
                x1, x2 = a, b
            x1 *= abs(scale)
            x2 *= abs(scale)
            print(x1, x2)
            if self.xy.index(index) == 0:
                self.ax.set_xlim(x1, x2)
            else:
                self.ax.set_ylim(x1, x2)

    def set_xlim(self, x1=None, x2=None):
        self._set_lim(x1, x2, 0)

    def set_ylim(self, y1=None, y2=None):
        self._set_lim(y1, y2, 1)

    def set_zlim(self, z1=None, z2=None):
        self._set_lim(z1, z2, 2)

    def set_lim(self, lim):
        self.set_xlim(-lim, lim)
        self.set_ylim(-lim, lim)
        self.set_zlim(-lim, lim)

    def get_xy(self, p):
        x = p[self.xy[0]] * np.abs(self.scale[0])
        y = p[self.xy[1]] * np.abs(self.scale[1])
        return x, y

    def _draw(self, *primitives):
        for pp in primitives:
            childs = pp.render()
            if childs is pp:
                getattr(self, "draw_" + pp.__class__.__name__)(pp)
            else:
                self.draw(*childs)

    def draw(self, *primitives):
        for pp in primitives:
            self.objects.add(pp)
        return self._draw(*primitives)

    def redraw(self):
        self.ax.clear()
        self._draw(*self.objects)

    def draw_Point(self, point):
        x, y = self.get_xy(point.p)
        self.ax.plot([x], [y], color=point.color, marker=point.marker)

    def draw_Line(self, line):
        x1, y1 = self.get_xy(line.p1)
        x2, y2 = self.get_xy(line.p2)
        self.ax.plot(
            [x1, x2],
            [y1, y2],
            color=line.color,
            marker=line.marker,
            linestyle=line.linestyle,
        )

    def draw_Text(self, text):
        x, y = self.get_xy(text.center)
        self.ax.text(x, y, text.text)

    def draw_Curve(self, curve):
        x, y = self.get_xy(curve.xyz)
        self.ax.plot(
            x, y, color=curve.color, marker=curve.marker, linestyle=curve.linestyle
        )


def rot_s(psi):
    """Rotation around z axis. Positive angle move x into y."""
    c = np.cos(psi)
    s = np.sin(psi)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rot_x(phi):
    """Rotation around x axis. Positive angle move y into z."""
    c = np.cos(phi)
    s = np.sin(phi)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rot_y(theta):
    """Rotation around y axis. Positive angle move z into x."""
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def rot_mad(theta=0, phi=0, psi=0):
    """Intrinsic rotaion. Order psi(roll), phi(pitch), theta(jaw)"""
    return rot_y(theta) @ rot_x(-phi) @ rot_s(psi)


def get_w_from_angles(theta, phi, psi, w):
    costhe = np.cos(theta)
    cosphi = np.cos(phi)
    cospsi = np.cos(psi)
    sinthe = np.sin(theta)
    sinphi = np.sin(phi)
    sinpsi = np.sin(psi)
    w[0, 0] = +costhe * cospsi - sinthe * sinphi * sinpsi
    w[0, 1] = -costhe * sinpsi - sinthe * sinphi * cospsi
    w[0, 2] = sinthe * cosphi
    w[1, 0] = cosphi * sinpsi
    w[1, 1] = cosphi * cospsi
    w[1, 2] = sinphi
    w[2, 0] = -sinthe * cospsi - costhe * sinphi * sinpsi
    w[2, 1] = +sinthe * sinpsi - costhe * sinphi * cospsi
    w[2, 2] = costhe * cosphi
    return w


def get_global_orbit(xs, ys, zs, theta, phi, psi, x, y, s):
    w = get_w_from_angles(theta, phi, psi, np.zeros([3, 3]))
    return np.array([xs, ys, ss]) + w @ np.array([x, y, s])


def advance(v, w, length=0, angle=0, tilt=0):
    if angle == 0:
        dz = np.array([0, 0, length])
        return v + w @ dz, w
    else:
        rho = length / angle
        ca = np.cos(angle)
        sa = np.sin(angle)
        ct = np.cos(tilt)
        st = np.sin(tilt)
        vv = np.array([rho * (ca - 1), 0, rho * sa])
        rot_y = np.array([[ca, 0, -sa], [0, 1, 0], [sa, 0, ca]])
        rot_z = np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]])
        roti_z = np.array([[ct, st, 0], [-st, ct, 0], [0, 0, 1]])
        vv = rot_z @ vv
        ss = rot_z @ rot_y @ roti_z
        return np.dot(w, vv) + v, np.dot(w, ss)


class Reference:
    @classmethod
    def from_mad(cls, x=0, y=0, z=0, theta=0, phi=0, psi=0):
        v = np.array([x, y, z])
        w = rot_mad(theta, phi, psi)
        return cls(v, w)

    def __init__(self, v, w):
        self.v = np.array(v, dtype=float)
        self.w = np.array(w, dtype=float)

    def advance(self, length=0, angle=0, tilt=0):
        v, w = advance(self.v, self.w, length=length, angle=angle, tilt=tilt)
        return Reference(v, w)

    def transform(self, p):
        "transform local p in global coordinates (x,y,s)"
        return self.v + self.w @ np.array(p)

    def curve_ahead(self, length=0, angle=0, tilt=0, ds=0.01):
        steps = int(length / ds)
        vv = np.zeros((3, steps + 1))
        ds = length / steps
        da = angle / steps
        v = self.v
        w = self.w
        vv[:, 0] = v
        for ss in range(steps):
            v, w = advance(v, w, ds, da, tilt)
            vv[:, ss + 1] = v
        return vv

    def render(self):
        o = Point(self.transform((0, 0, 0)), color="k")
        x = Point(self.transform((1, 0, 0)), color="r")
        y = Point(self.transform((0, 1, 0)), color="g")
        z = Point(self.transform((0, 0, 1)), color="b")
        xv = Line(o, x, color="r")
        yv = Line(o, y, color="g")
        zv = Line(o, z, color="b")
        lx = Text("X", x)
        ly = Text("Y", y)
        lz = Text("S", z)
        return [o, x, y, z, xv, yv, zv, lx, ly, lz]


class Survey:
    mad_columns = [
        "name",
        "s",
        "l",
        "angle",
        "x",
        "y",
        "z",
        "theta",
        "phi",
        "psi",
        "globaltilt",
        "mech_sep",
        "assembly_id",
        "slot_id",
    ]

    @classmethod
    def from_data(cls, data):
        self = cls()
        for cc in cls.mad_columns:
            setattr(self, cc, data[cc])
        return self

    @classmethod
    def from_mad_table(cls, mad):
        return cls.from_data(mad.table.survey)

    def get_global_orbit(self, idx, x, y, s):
        vv = p.array([self.x[idx], self.y[idx], self.z[idx]])
        ww = np.zeros((3, 3))
        get_w_from_angles(self.theta[idx], self.phi[idx], self.psi[idx], ww)
        ll = np.zeros([x, y, s])
        return vv + ww @ ll

    @classmethod
    def open(cls, filename):
        return cls.from_data(tfsdata.open(filename))


class MadElem:
    def __init__(self, name, length=0, angle=0, tilt=0, ref=None):
        self.name = name
        self.length = length
        self.angle = angle
        self.tilt = tilt
        if ref is None:
            ref = Reference.from_mad()
        self.set_ref(ref)

    def set_ref(self, start):
        self.ref_start = start
        self.ref_middle = start.advance(
            length=self.length / 2, angle=self.angle / 2, tilt=self.tilt
        )
        self.ref_end = start.advance(
            length=self.length, angle=self.angle, tilt=self.tilt
        )
        return self

    def render(self, ends=True):
        ps = Point(self.ref_start.v, color="b")
        pm = Point(self.ref_middle.v, color="r")
        pe = Point(self.ref_end.v, color="b")
        cav = self.ref_start.curve_ahead(
            self.length / 2, angle=self.angle / 2, tilt=self.tilt
        )
        ca = Curve(cav, color="k")
        cbv = self.ref_middle.curve_ahead(
            self.length / 2, angle=self.angle / 2, tilt=self.tilt
        )
        cb = Curve(cbv, color="k")
        tm = Text(self.name, pm, color="k")
        return [ca, cb, ps, pm, pe, tm]


class MadLine:
    def __init__(self, *elements, ref=None):
        if ref is None:
            ref = Reference.from_mad()
        self.ref_start = ref
        self.elements = list(elements)

    def add(self, name, length=0, angle=0, tilt=0):
        if len(self.elements) == 0:
            ref = self.ref_start
        else:
            ref = self.elements[-1].ref_end
        el = MadElem(name, length=length, angle=angle, tilt=tilt, ref=ref)
        self.elements.append(el)
        return el

    def render(self):
        return self.elements
