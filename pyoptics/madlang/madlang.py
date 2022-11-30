from .madobj import Elem, Expr, ExprList, Sequence, Line

from .parser import parse, pyname, parses


def fromast(ast, name="mad", root=None, special=None):
    if special is None:
        special = []
    if root is None:
        root = Elem(name=name, _orig={}, _rorig={})
    macros = {}
    current_seq = None
    for st in ast:
        if st is None:
            continue
        # print 'st',st
        madname, value = st
        name = pyname(madname)
        kind = value[0]
        if hasattr(root, "_orig"):
            root._orig[madname] = name
        if hasattr(root, "_rorig"):
            root._rorig[name] = madname
        if kind == "variable":
            # print 'value',value
            if name in special:
                val = evaluate(value[1], None)
            else:
                val = evaluate(value[1], root)
            # print 'val',val
            root[name] = val
        elif kind == "expression":
            value = mkexpr(value[1])
            root[name] = value
        elif kind == "element":
            proto = value[1][0]
            if proto is not None:
                proto = pyname(proto)
            if proto == "macro":
                out = []
                for st in ast:
                    out.append(st)
                    if st == ("statement", ("value", "};")):
                        break
                macros[name] = out
            elif name == "endsequence":
                current_seq = None
            elif name == "return":
                break
            elif proto == "sequence":
                ne = Sequence(madname, parent=root["sequence"])
                current_seq = ne
            elif proto == "line":
                # right hand side of line statement is not interpret as assignment
                # later by fromast(attrs,....)
                ne = Line(name=madname, parent=root["line"], value=value[1][1][1])
            elif proto == None:
                ne = root[name]
            else:
                proelem = root[proto]
                ne = Elem(name=madname, parent=proelem)
            root[name] = ne
            attrs = value[1][1:]
            fromast(attrs, name=madname, root=ne, special=["refer", "From", "apertype"])
            if current_seq is not None and current_seq is not ne:
                current_seq.append(name, ne)
    # print macros
    return root


def load(fh, name=None, root=None):
    return fromast(parse(fh), name=name, root=root)


def loads(s, root=None):
    return fromast(parses(s), root=root)


_pyopen = open


def open(fn, root=None):
    if fn.endswith(".gz"):
        fh = gzip.open(fn)
    else:
        fh = _pyopen(fn)
    if root is not None:
        fn = "%s,%s" % (root.name, fn)
    return load(fh, name=fn, root=root)


def evaluate(value, lcl):
    if type(value) is str:
        try:
            value = pyname(value)
            if lcl:
                value = eval(value, Elem.gbl, lcl)
        except NameError as e:
            print("Warning", value, "not evaluated")
            print(e)
    elif type(value) is list:
        value = [evaluate(i, lcl) for i in value]
    return value


def mkexpr(value):
    if type(value) is str:
        value = pyname(value)
        value = Expr(value)
    elif type(value) is list:
        value = [pyname(v) for v in value]
        value = ExprList(*value)
    return value
