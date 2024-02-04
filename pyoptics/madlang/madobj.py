from collections import namedtuple
import numpy as np

import ast


def get_names(expr):
    tree = ast.parse(expr)
    out = []
    toremove = []
    for node in ast.walk(tree):
        if type(node) is ast.Attribute:
            name = _get_names(node.value).pop()
            out.append(".".join([name, node.attr]))
            toremove.append(name)
        elif type(node) is ast.Name:
            out.append(node.id)
    for name in toremove:
        out.remove(name)
    return set(out)


class ExprUndefined(Exception):
    pass


def get_attrs(obj):
    import types

    if not hasattr(obj, "__dict__"):
        return []  # slots only
    # if not isinstance(obj.__dict__, (dict, types.DictProxyType)):
    if not isinstance(obj.__dict__, dict):
        raise TypeError("%s.__dict__ is not a dictionary" "" % type(obj.__name__))
    return list(obj.__dict__.keys())


def dir2(obj):
    attrs = set()
    if not hasattr(obj, "__bases__"):
        # obj is an instance
        if not hasattr(obj, "__class__"):
            # slots
            return sorted(get_attrs(obj))
        klass = obj.__class__
        attrs.update(get_attrs(klass))
    else:
        # obj is a class
        klass = obj
    for cls in klass.__bases__:
        attrs.update(get_attrs(cls))
        attrs.update(dir2(cls))
    attrs.update(get_attrs(obj))
    return list(attrs)


class Elem(object):
    #  __slots__=['name','parent','_data','_ns','_orig']
    gbl = {}

    def __init__(self, name=None, parent=None, _orig=None, _rorig=None, **kwargs):
        self._data = {}
        self._ns = Elem.gbl
        self._parent = None
        if parent is not None:
            self._parent = parent.name
            self._data.update(parent._data)
        self.name = name
        if _orig is not None:
            self._orig = _orig
        if _rorig is not None:
            self._rorig = _rorig
        self._data.update(kwargs)

    def __repr__(self):
        def _str(k, v):
            if hasattr(v, "_get"):
                vv = "%s -> %s" % (v, self[k])
            else:
                vv = "%s" % (v)
            return "%-15s: %s" % (k, vv)

        if len(self._data) < 60:
            data = "\n%s\n" % "\n".join(
                _str(k, v) for k, v in sorted(self._data.items())
            )
        else:
            data = "%d objects" % len(self._data)
        if self._parent is not None:
            data = "%s,%s" % (self._parent, data)
        tmp = "<%s: %s>" % (self.name, data)
        return tmp

    def __getitem__(self, k, opt=None):
        try:
            v = self._data[k]
        except KeyError:
            if self._ns is not None:
                v = self._ns[k]
            else:
                raise KeyError("%s not found" % (k,))
        if hasattr(v, "_get"):
            v = v._get(self, self.gbl)
        return v

    def __setitem__(self, k, v):
        self._data[k] = v
        if hasattr(v, "_bind"):
            v._bind(self, k)

    def __delitem__(self, k):
        v = self[k]
        del self._data[k]

    def _bind(self, obj, k):
        self._ns = obj

    def __contains__(self, k):
        return k in self._data

    def __getattribute__(self, k):
        if k.startswith("_"):
            return object.__getattribute__(self, k)
        else:
            try:
                return self[k]
            except KeyError:
                return object.__getattribute__(self, k)

    #    print k, k in self.__class__.__dict__
    #    if k in self.__class__.__dict__ or k.startswith('_'):
    #       return object.__getattribute__(self,k)
    #    elif k in self._data:
    #        return self[k]
    #    else:
    #      raise AttributeError('%s missing in %s'%(k,self))

    def __setattr__(self, k, v):
        if k.startswith("_"):
            self.__dict__[k] = v
        else:
            self._data[k] = v

    def __delattr__(self, k):
        if k.startswith("_"):
            del self.__dict__[k]
        else:
            del self._data[k]

    def __dir__(self):
        #        out = dir2(self)
        out = list(self.__dict__.keys())
        out.extend(list(self._data.keys()))
        return out

    def __iter__(self):
        return self._data.__iter__()

    def __hasattr__(self, k):
        return k in self._data

    def build_dep(self):
        out = {}
        for key, att in list(self._data.items()):
            if hasattr(att, "expr"):
                for name, idx, expr in att.get_names():
                    out.setdefault(name, []).append((key, idx, expr))
            elif hasattr(att, "_data"):
                for name, lst in list(att.build_dep().items()):
                    for ll in lst:
                        out.setdefault(name, []).append((key,) + ll)
        return out

    def print_dep(self, key, indent="", deps=None):
        if deps is None:
            deps = self.build_dep()
        for dep in deps.get(key, []):
            expr = dep[-1]
            idx = dep[-2]
            objs = ".".join(dep[:-2])
            if idx is not None:
                objs = "%s[%d]" % (objs, idx)
            print("%s%-20s = %s" % (indent, objs, expr))
            if objs in deps:
                self.print_dep(objs, indent + "  ", deps)

    def __call__(self, expr):
        return eval(expr, {}, self)

    def vars_to_dict(self):
        out = {}
        for k, v in self._data.items():
            if np.isscalar(v):
                out[self._rorig[k]] = v
        return out


class Expr(object):
    __slots__ = "expr"

    def __init__(self, expr):
        self.expr = expr
        if expr.endswith(";"):
            self.expr = compile(expr, expr, "exec")
        else:
            self.expr = compile(expr, expr, "eval")

    def _get(self, lcl, gbl={}):
        try:
            return eval(self.expr, gbl, lcl)
        except NameError as e:
            print(("Warning %r undefined" % (self)))
            # return ExprUndefined(e.message)
            return 0

    def __repr__(self):
        return "Expr(%r)" % (self.expr.co_filename)

    def get_names(self, ix=None):
        names = self.expr.co_names
        expr = [self.expr.co_filename] * len(names)
        return list(zip(names, [ix] * len(names), expr))


class ExprList(object):
    __slots__ = "expr"

    def __init__(self, *expr):
        self.expr = [Expr(ex) for ex in expr]

    def _get(self, lcl, gbl={}):
        return [ex._get(lcl, gbl) for ex in self.expr]

    def __repr__(self):
        exp = ",".join(["%r" % ex.expr.co_filename for ex in self.expr])
        return "ExprList(%s)" % exp

    def get_names(self):
        out = []
        for ix, ex in enumerate(self.expr):
            out.extend(ex.get_names(ix))
        return out


class SeqElem(namedtuple("SeqElem", "at From mech_sep slot_id")):
    @classmethod
    def from_dict(cls, data):
        return cls(*[data[k] for k in cls._fields])


classes = dict(
    Drift=namedtuple("Drift", "length"),
    DriftExact=namedtuple("DriftExact", "length"),
    Multipole=namedtuple("Multipole", "knl ksl hxl hyl length"),
    Cavity=namedtuple("Cavity", "voltage frequency lag"),
    XYShift=namedtuple("XYShift", "dx dy"),
    SRotation=namedtuple("SRotation", "angle"),
    Line=namedtuple("Line", "elems"),
    BeamBeam4D=namedtuple(
        "BeamBeam4D",
        " ".join(
            [
                "q_part",
                "N_part",
                "sigma_x",
                "sigma_y",
                "beta_s",
                "min_sigma_diff",
                "Delta_x",
                "Delta_y",
            ]
        ),
    ),
    # BeamBeam6D=namedtuple('BeamBeam6D', ' '.join([
    #         'q_part N_part_tot sigmaz N_slices min_sigma_diff threshold_singular',
    #         'phi alpha',
    #         'Sig_11_0 Sig_12_0 Sig_13_0',
    #         'Sig_14_0 Sig_22_0 Sig_23_0',
    #         'Sig_24_0 Sig_33_0 Sig_34_0 Sig_44_0'
    #         'delta_x delta_y'
    #         'x_CO px_C0 y_CO py_CO sigma_CO delta_CO'
    #         'Dx_sub Dpx_sub Dy_sub Dpy_sub Dsigma_sub Ddelta_sub'
    #         'enabled']))
    BeamBeam6D=namedtuple("BeamBeam6D", "BB6D_data"),
)


def flatten(lst):
    if isinstance(lst, (list, tuple)):
        for i in lst:
            if isinstance(i, (list, tuple)):
                for j in flatten(i):
                    yield j
            else:
                yield i
    else:
        yield lst


class fakedict(object):
    def __getitem__(self, k):
        return k


fakedict = fakedict()


class Line(Elem):
    def __init__(self, name=None, parent=None, value=None, **kwargs):
        Elem.__init__(self)
        self.value = value

    def flatten_names(self):
        return tuple(flatten(eval(self.value, {}, fakedict)))

    def flatten_objects(self):
        def _flatten(lst):
            if isinstance(lst, (list, tuple)):
                for elem in lst:
                    if type(elem) is tuple:
                        for ee in _flatten(elem):
                            yield ee
                    elif elem.keyword == "line":
                        for ee in elem.flatten_objects():
                            yield ee
                    else:
                        yield elem
            else:
                yield lst

        elems = list(_flatten(eval(self.value, {}, self._ns)))
        return elems

    def expand_struct(self, convert=classes):
        Drift = convert["Drift"]
        DriftExact = convert["DriftExact"]
        Multipole = convert["Multipole"]
        Cavity = convert["Cavity"]
        # XYshift = convert['XYShift']
        # SRotation = convert['SRotation']
        rest = []
        names = self.flatten_names()
        elems = self.flatten_objects()
        newelems = []
        types = []
        iconv = []
        icount = 0
        for elem in elems:
            if elem.keyword == "drift":
                newelems.append(DriftExact(length=elem.l))
                types.append("DriftExact")
                icount += 1
            elif elem.keyword == "Multipole":
                newelems.append(
                    Multipole(
                        knl=elem.knl,
                        ksl=elem.ksl,
                        length=elem.lrad,
                        hxl=elem.knl[0],
                        hyl=elem.ksl[0],
                    )
                )
                types.append(elem.keyword)
                icount += 1
            elif elem.keyword in ["hkicker"]:
                ne = Multipole(
                    knl=[-elem.kick], ksl=[], length=elem.lrad, hxl=elem.kick, hyl=0
                )
                newelems.append(ne)
                types.append("Multipole")
                icount += 1
            elif elem.keyword in ["vkicker"]:
                ne = Multipole(
                    knl=[], ksl=[elem.kick], length=elem.lrad, hxl=0, hyl=elem.kick
                )
                newelems.append(ne)
                types.append("Multipole")
                icount += 1
            elif elem.keyword in ["rfcavity"]:
                nvolt = elem.volt * 1e6
                ne = Cavity(
                    voltage=nvolt, frequency=elem.freq * 1e6, lag=elem.lag * 360
                )
                newelems.append(ne)
                types.append("Cavity")
                icount += 1
            else:
                rest.append(el)
            iconb.append(icount)
        newelems = [dict(i._asdict()) for i in newelems]
        return list(zip(names, types, newelems)), rest, iconv


class Sequence(Elem):
    _fields = "at From mech_sep slot_id".split()

    def __init__(self, name=None, parent=None, elems=None, **kwargs):
        Elem.__init__(self)
        if elems is None:
            self.elems = []

    def get_names(self):
        return [el.name for el in self.elems]

    def append(self, name, elem):
        ne = Elem(
            name=name,
            at=elem._data.get("at"),
            From=elem._data.get("From"),
            slot_id=elem._data.get("slot_id"),
            mech_sep=elem._data.get("mech_sep"),
        )
        self.elems.append(ne)

    def expand_struct(self, convert=classes):
        elems = []
        rest = []
        count = {}
        Drift = convert["Drift"]
        DriftExact = convert["DriftExact"]
        Multipole = convert["Multipole"]
        Cavity = convert["Cavity"]
        # XYshift = convert['XYShift']
        # SRotation = convert['SRotation']
        pos = {}
        lasts = 0
        drifts = {}
        for el in self.elems:
            if el.From is None:
                pos[el.name] = el.at
        iconv = []
        icount = 0
        for el in self.elems:
            elem = self[el.name]
            s = el.at
            if el.From is not None:
                s += pos[el.From]
            l = elem.l
            ldrift = s - l / 2 - lasts
            if lasts + l / 2 < s:
                dr = drifts.get(ldrift)
                if dr is None:
                    drname = "drift_%d" % len(drifts)
                    dr = DriftExact(length=ldrift)
                    drifts[l] = dr
                elems.append((drname, "DriftExact", dr))
                icount += 1
            if elem.keyword == "multipole":
                ne = Multipole(
                    knl=elem.knl,
                    ksl=elem.ksl,
                    length=elem.lrad,
                    hxl=elem.knl[0],
                    hyl=elem.ksl[0],
                )
                elems.append((elem.name, "Multipole", ne))
                icount += 1
            elif elem.keyword in [
                "marker",
                "hmonitor",
                "vmonitor",
                "instrument",
                "monitor",
                "rcollimator",
            ]:
                ne = Drift(length=0)
                elems.append((elem.name, "Drift", ne))
                icount += 1
            elif elem.keyword in ["hkicker"]:
                ne = Multipole(
                    knl=[-elem.kick], ksl=[], length=elem.lrad, hxl=elem.kick, hyl=0
                )
                elems.append((elem.name, "Multipole", ne))
                icount += 1
            elif elem.keyword in ["vkicker"]:
                ne = Multipole(
                    knl=[], ksl=[elem.kick], length=elem.lrad, hxl=0, hyl=elem.kick
                )
                elems.append((elem.name, "Multipole", ne))
                icount += 1
            elif elem.keyword in ["rfcavity"]:
                nvolt = elem.volt * 1e6
                ne = Cavity(
                    voltage=nvolt, frequency=elem.freq * 1e6, lag=elem.lag * 360
                )
                elems.append((elem.name, "Cavity", ne))
                icount += 1
            else:
                rest.append((elem.name, elem))
            iconv.append(icount)
            lasts = s + l / 2
        return elems, rest, iconv


elements = [
    "hkicker",
    "vkicker",
    "ecollimator",
    "instrument",
    "kicker",
    "marker",
    "monitor",
    "multipole",
    "octupole",
    "quadrupole",
    "rbend",
    "sbend",
    "rcollimator",
    "rfcavity",
    "sextupole",
    "solenoid",
    "tkicker",
    "placeholder",
    "drift",
    "hmonitor",
    "vmonitor",
    "sequence",
    "collimator",
    "dipedge",
]

commands = [
    "beam",
    "endmatch",
    "global",
    "jacobian",
    "match",
    "track",
    "twiss",
    "use",
    "value",
    "vary",
    "line",
    "track",
    "start",
    "run",
    "endtrack",
    "set",
]


for e in elements:
    Elem.gbl[e] = Elem(e, keyword=e, l=0)


Elem.gbl["multipole"]._data.update(knl=[0], ksl=[0], lrad=0)
Elem.gbl["rfcavity"]._data.update(lag=0)
Elem.gbl["hkicker"]._data.update(lrad=0)
Elem.gbl["vkicker"]._data.update(lrad=0)
Elem.gbl["kicker"]._data.update(lrad=0)

for e in commands:
    Elem.gbl[e] = Elem(e)

import math

for k in dir(math):
    if not k.startswith("_"):
        Elem.gbl[k] = getattr(math, k)

Elem.gbl["sinc"] = lambda x: math.sin(math.pi * x) / (math.pi * x)
Elem.gbl["twopi"] = 2 * math.pi

names = ["proton", "true", "false", "electron", "ion"]

for nn in names:
    Elem.gbl[nn] = nn
