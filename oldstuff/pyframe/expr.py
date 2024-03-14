class expr(object):
    __slots__ = ["_doc", "expr", "code", "lcl"]

    def __init__(self, expr, lcl=None, _doc=None):
        self.lcl = lcl
        self._doc = _doc
        self._update(expr)

    def _call(self, *argl, **argd):
        return eval(self.code, globals(), self.lcl)

    def _update(self, expr, *args):
        self.expr = expr
        if expr.endswith(";"):
            self.code = compile(expr, expr, "exec")
        else:
            self.code = compile(expr, expr, "eval")

    def _env(self, lcl, *args):
        if self.lcl is None:
            self.lcl = lcl

    def __repr__(self):
        return "%s -> %s" % (self.expr, self._call())


if __name__ == "__main__":
    import doctest

    doctest.testmod()
