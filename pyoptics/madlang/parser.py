import re

comment = [
    re.compile(r"//.*"),
    re.compile("/\*.*?\*/", re.S),
    re.compile(r"!.*"),
    re.compile(r"real |const |shared ", re.I),
    re.compile(r"\s", re.S),
]
statement = re.compile("[^;]*;")
variable = re.compile("(^[\w\.]+)(:?=)([^;:,]+)")
element = re.compile("(^[\w\.()]+)(:([\w\.]+))?([,=]?.+)?")
ropt = re.compile("(^[\w]+)")


def parse_attr(attr):
    ch = False
    out = []
    if attr:
        for i in attr:
            if i == "{" or i == "(":
                ch = True
            if i == "}" or i == ")":
                ch = False
            if ch and (i == ","):
                out.append(" ")
            elif i != ";":
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
        name, proto, attrlist = res[0], res[2], res[3]
        value = []
        res = ["element", value]
        value.append(proto)
        if attrlist:
            for i in parse_attr(attrlist):
                attr = parse_variable(i)
                if attr is None:
                    attr = [
                        "value",
                        i.replace(" ", ","),
                    ]  # use trick in parseattr ',' -> ' '
                value.append(attr)
        return [name, res]


def parses(s):
    s = s.lower()
    for r in comment:
        s = r.sub("", s)
    s = s.lower()
    current_sequence = None
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
            res = ("statement", ("value", st))
        yield res


def parse(fh):
    return parses(fh.read())


def no_dots(x):
    return x.group().replace(".", "_")


madname = re.compile(r"([a-z_][a-z_0-9\.]*)")


def pyname(n):
    n = n.lower()
    n = madname.sub(no_dots, n)
    n = n.replace("^", "**")
    n = n.replace("->", ".")
    if n == "from":
        n = "From"
    return n
