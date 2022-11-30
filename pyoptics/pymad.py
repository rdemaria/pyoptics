from collections import namedtuple, defaultdict
import re

from numpy import *

from .madxseqdata import madxseqdata


gbl = globals()


def quad(l, kn0l, ks0l, kn1l, ks1l, br=1, gr=7000):
    r = zeros((6, 6))
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
            skx = sqrt(kx)
            cx = cos(skx * l)
            sx = sin(skx * l) / skx
            dx = (1 - cx) / kx
            fx = (l - sx) / kx
        else:
            skx = sqrt(-kx)
            cx = cosh(skx * l)
            sx = sinh(skx * l) / skx
            dx = (1 - cx) / kx
            fx = (l - sx) / kx
        ky = -kn1l / l + (ks0l / l) ** 2
        if abs(ky) < 1e-10:
            cy = 1 - l**2 * ky / 2
            sy = l - l**3 * ky / 6
            dy = l**2 / 2
        elif ky > 0:
            sky = sqrt(ky)
            cy = cos(sky * l)
            sy = sin(sky * l) / sky
            dy = (1 - cy) / ky
        else:
            sky = sqrt(-ky)
            cy = cosh(sky * l)
            sy = sinh(sky * l) / sky
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


def propbeta(r, betx, alfx, mux, bety, alfy, muy):
    t1 = r[0, 0] * betx - r[0, 1] * alfx
    t2 = r[1, 0] * betx - r[1, 1] * alfx
    newbetx = (t1**2 + r[0, 1] ** 2) / betx
    newalfx = -(t1 * t2 + r[0, 1] * r[1, 1]) / betx
    newmux = mux + arctan2(r[0, 1], t1) / (2 * pi)
    t3 = r[2, 2] * bety - r[2, 3] * alfy
    t4 = r[3, 2] * bety - r[3, 3] * alfy
    newbety = (t3**2 + r[2, 3] ** 2) / bety
    newalfy = -(t3 * t4 + r[2, 3] * r[3, 3]) / bety
    newmuy = muy + arctan2(r[2, 3], t3) / (2 * pi)
    return [newbetx, newalfx, newmux, newbety, newalfy, newmuy]


class expr(object):
    def __init__(self, expr):
        self.expr = expr
        self.code = compile(expr, "eval", "eval")

    def _on_get(self, dct, name):
        return eval(self.code, gbl, dct)


class Frame(object):
    _parent = {}
    _names = {}

    def __init__(self, **names):
        self.__dict__["_names"] = {}
        self._names.update(names)

    def __getitem__(self, k):
        if k in self._names:
            v = self._names[k]
        else:
            v = self._parent[k]
        if hasattr(v, "_on_get"):
            return v._on_get(self, k)
        else:
            return v

    def __setitem__(self, k, v):
        if hasattr(v, "_on_set"):
            v._on_set(self, k)
        self._names[k] = v

    def __delitem__(self, k):
        del self._names[k]

    def keys(self):
        return list(self._names.keys())

    def copy(self):
        return self.__class__(self._names.copy())

    def _on_set(self, setter, k):
        self.__dict__["_parent"] = setter

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError

    def __setattr__(self, k, v):
        self[k] = v

    def __delattr__(self, k):
        del self[k]


class Elem(Frame):
    def __init__(self, l=0, kn0l=0, ks0l=0, kn1l=0, ks1l=0, **names):
        names = dict(l=l, kn0l=kn0l, ks0l=ks0l, kn1l=kn1l, ks1l=ks1l, **names)
        self.__dict__["_names"] = names

    def rmatrix(self, br=1, gr=1):
        return quad(self.l, self.kn0l, self.ks0l, self.kn1l, self.ks1l, br, gr)

    def track(self, init):
        r = self.rmatrix(init.br, init.gr)
        s = init.s + self.l
        res = propbeta(
            r, init.betx, init.alfx, init.mux, init.bety, init.alfy, init.muy
        )
        betx, alfx, mux, bety, alfy, muy = res
        disp = array([[init.dx, init.dpx, init.dy, init.dpy, 0, 1]]).T
        newdisp = dot(r, disp)
        out = [s] + res + list(dot(r, disp)[0:4, 0])
        return Beam(*out)


class Line(Frame):
    def __init__(self, *elems):
        self.__dict__["elems"] = list(elems)
        self.__dict__["_names"] = {}

    def track(self, init):
        data = {"start": init}
        for name in self.elems:
            el = self[name]
            init = el.track(init)
            data[name] = init
        return TwissTable(["start"] + self.elems, data)


class Beam(Frame):
    _default = dict(
        s=0,
        betx=1,
        alfx=0,
        mux=0,
        bety=1,
        alfy=0,
        muy=0,
        dx=0,
        dpx=0,
        dy=0,
        dpy=0,
        br=1,
        gr=7000,
    )
    cols = "s betx alfx mux bety alfy muy dx dpx dy dpy br gr".split()

    def __init__(self, *vals, **names):
        _names = Beam._default.copy()
        for k, v in zip(Beam.cols, vals):
            _names[k] = v
        _names.update(**names)
        self.__dict__["_names"] = _names


class TwissTable(Frame):
    _fmt_str = "%-8s"
    _fmt_dbl = "%8.3f"

    def __init__(self, seq, data):
        self.__dict__["_names"] = data
        self.__dict__["_seq"] = seq

    def print_table(self, cols=Beam.cols):
        cols = cols.split()
        print(" ".join(["%-8s" % n for n in cols]))
        for name in self._seq:
            line = []
            for col in cols:
                if col == "name":
                    line.append(self._fmt_str % name)
                else:
                    line.append(self._fmt_dbl % self[name][col])
            print(" ".join(line))


twopi = 2 * pi

comment = [
    re.compile(r"//.*"),
    re.compile("/\*.*?\*/", re.S),
    re.compile(r"!.*"),
    re.compile(r"real |const |shared ", re.I),
    re.compile(r"\s", re.S),
]
statement = re.compile("[^;]*;")
variable = re.compile("(^[\w\.]+)(:?=)([^;:,]+)")
element = re.compile("(^[\w\.]+)(:([\w\.]+))?(,.+)?")
ropt = re.compile("(^[\w]+)")


def parse_attr(attr):
    ch = False
    out = []
    if attr:
        for i in attr:
            if i == "{":
                ch = True
            if i == "}":
                ch = False
            if ch and (i == ","):
                out.append(" ")
            else:
                out.append(i)
    return out and ("".join(out[1:])).split(",") or out


def parse_variable(st):
    res = variable.match(st)
    if res:
        res = res.groups()
        cls = res[1] == ":=" and "expression" or "variable"
        name = res[0]
        expr = res[2]
        if expr.startswith("{"):
            expr = expr[1:-1].split()  # use trick in parseattr ',' -> ' '
        return [name, [cls, expr]]


def parse_element(st):
    res = element.match(st)
    if res:
        res = res.groups()
        name, proto, attr = res[0], res[2], res[3]
        value = []
        res = ["element", value]
        value.append(["proto", proto])
        if attr:
            for i in parse_attr(attr):
                value.append(parse_variable(i))
        return [name, res]


def parses(s):
    for r in comment:
        s = r.sub("", s)
    s = s.lower()
    current_Sequence = None
    for st in statement.finditer(s):
        st = st.group()
        while 1:
            res = parse_variable(st)
            if res:
                break
            res = parse_element(st)
            if res:
                break
            break
        if res is None:
            res = (("class", "statement"), ("value", st))
        yield res


def parse(fh):
    return parses(fh.read())


basenames = [
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
]


base = Frame()
for i in basenames:
    base[i] = Frame(name=i, elemclass=i)


def fromast(ast, root=None, lcl=None, special=[], name="mad"):
    if root is None:
        root = Frame(name=name)
    if lcl is None:
        lcl = root
    current_seq = None
    root._originalnames = {}
    for i, st in enumerate(ast):
        if st is None:
            continue
        name, value = st
        nname = pyname(name)
        root._originalnames[nname] = name
        name = nname
        kind = value[0]
        if kind == "variable":
            if name in special:
                value = evaluate(value[1], None)
            else:
                value = evaluate(value[1], lcl)
            setattr(root, name, value)
        elif kind == "expression":
            value = mkexpr(value[1])
            setattr(root, name, value)
        elif kind == "element":
            # [element, [[proto,.],[.,.]]] -> [[proto,.],[.,..]]]
            value = value[1]
            proto = value[0][1]  # can be None
            if proto:
                proto = pyname(proto)
            if name == "endSequence":
                current_seq.n_elems = len(current_seq._elements)
                current_seq = None
                current_pos = 0
            elif name == "return":
                pass
            else:
                if name in root:  # element already defined
                    ne = getattr(root, name)
                elif proto == "Sequence":  # special element sequence
                    ne = Sequence(name)
                    ne._initseq()
                    current_seq = ne
                else:  # element not defined, therefore has proto
                    proto = getattr(root, proto)
                    ne = Frame(name, proto)
                setattr(root, name, ne)
                ne._parent = root
                fromast(
                    value[1:], root=ne, lcl=lcl, special=["refer", "From", "apertype"]
                )
                if current_seq is not None and current_seq is not ne:
                    current_seq._append(ne)
    return root


class Sequence(Frame):
    refer = "centre"

    def _initseq(self):
        self._elements = []
        self._element_dict = defaultdict(list)

    def _append(self, newelem):
        if "From" in newelem:
            frm = newelem["From"]
        else:
            frm = None
        pos = len(self._elements)
        data = madxseqdata(pos, newelem, newelem.at, frm, self)
        self._elements.append(data)
        self._element_dict[newelem.name].append(data)

    def get_data(self, name):
        return [data for data in self._element_dict[name]]

    def get_pos(self, name):
        return [data.position for data in self.get_data(name)]

    def get_table(self, start=None, end=None, cols="name s l angle k1"):
        if start:
            start = self.get_data(start)[0].position
        else:
            start = 0
        if end:
            end = self.get_data(end)[-1].position + 1
        else:
            end = -1
        elems = self._elements[start:end]
        cols = cols.split()
        if "s" in cols:
            for i in elems:
                i.element.s = self.startpos(i)
        table = [[getattr(i.element, c, 0.0) for c in cols] for i in elems]
        table = list(zip(*table))
        table = list(map(array, table))
        res = namedtuple("result", cols)(*table)
        return res

    def get_elements(self, start=None, end=None):
        if start:
            start = self.get_data(start)[0].position
        else:
            start = 0
        if end:
            end = self.get_data(end)[-1].position + 1
        else:
            end = -1
        return [i.element for i in self._elements[start:end]]

    def startpos(self, data):
        if self.refer == "centre":
            position = data.at - data.element.l / 2.0
        else:
            raise self.refer + " Not implemented"
        if data.From:
            position += self.startpos(data.Sequence._element_dict[data.From][0])
        return position


def load(fh):
    return fromast(parse(fh))


def loads(fh):
    return fromast(parses(fh))


def no_dots(x):
    return x.group().replace(".", "_")


madname = re.compile(r"([a-z_][a-z_0-9\.]*)")


def pyname(n):
    n = n.lower()
    n = madname.sub(no_dots, n)
    n.replace("^", "**")
    if n == "from":
        n = "From"
    return n


def evaluate(value, lcl):
    if type(value) is str:
        try:
            value = pyname(value)
            if lcl:
                value = eval(value, globals(), lcl)
        except NameError:
            print("Warning", value, "not evaluated")
    elif type(value) is list:
        value = [evaluate(i, lcl) for i in value]
    return value


def mkexpr(value):
    if type(value) is str:
        value = pyname(value)
        value = expr(value)
    elif type(value) is list:
        value = mkexpr(repr(value))
    return value


# class Match(object):
#  def __init__(self,frame,varlist,fun):
#    pass
#  def __call__(self,x):
#    for name,val in zip(varlist,x):
#      frame[x]=val
#    return fun(frame)
#
# def ftosolve(f):
#  f.track(f.seq1)
#
# class Match(object):
#  def __init__(self,ftorun,frame,varlist,constraints):
#    self.ftorun=ftorun
#    self.varlist=varlist
#    self.constraints=constraints
#  def ftosolve(self,x):
#    bounds=[]
#    for (name,frame,ub,lb),val in zip(varlist,x):
#      frame[name]=val
#      bounds.append((ub,lb))
#    obj=self.ftorun(frame)
#    eq=[]
#    ineq=[]
#    for cons in self.constraints:
#      val=cons()
#      if cons.iseq():
#        eq.appen(val)
#      else:
#        ineq.appen(val)
#    return obj,eq,ineq,bounds
