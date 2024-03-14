from namedtuple import namedtuple

from elem import expr


class constraint(object):
    def __init__(self, expr, model, name="cons", weight=1, equality=1, active=True):
        self.name = name
        self.expr = expr
        self.model = model
        self.weight = weight
        self.equality = equality
        self.active = active
        self.mkexpr()

    def mkexpr(self):
        if "=" in self.expr:
            l = tuple(self.expr.split("="))
            self.expr = "%s-(%s)" % l
            self.equality = 1
        elif ">" in self.expr:
            l = tuple(self.expr.split(">"))
            self.expr = "-(%s)+(%s)" % l
            self.equality = 0
        elif "<" in self.expr:
            l = tuple(self.expr.split("<"))
            self.expr = "%s-(%s)" % l
            self.equality = 0
        self.code = expr(self.expr, lcl=self.model)


class vary(object):
    def __init__(
        self,
        model,
        pointer=None,
        name="cons",
        weight=1,
        min=None,
        max=None,
        slope=0,
        active=True,
    ):
        self.name = name
        self.pointer = pointer
        self.model = model
        self.weight = weight
        self.min = min
        self.max = max
        self.slope = slope
        self.getcons()

    def getcons(self):
        out = []
        if self.min:
            ex = "_self[%s]>%s" % (self.pointer, self.min)
            out.add(consntraint(self.model, ex, "%s_min" % self.name))
        if self.max:
            ex = "_self[%s]<%s" % (self.pointer, self.max)
            out.add(consntraint(self.model, ex, "%s_max" % self.name))


class match(object):
    vars = []  # name, weight, change
    cons = []  # name, code, weight, isboundary
    models = set()  # models to run in the macro

    def __getitem__(self, k):
        return self.__dict__[k]

    def __setitem__(self, k, v):
        self.__dict__[k] = v

    def __delitem__(self, k):
        del self.__dict__[k]

    class __metaclass__(type):
        def __setattr__(self, k, v):
            if isinstance(v, vary):
                v.name = k
                if v.pointer is None:
                    v.pointer = k
                self.vars.add(v)
                self.models.add(v.model)
            elif isinstance(v, constraint):
                v.name = k
                self.cons.add(v)
            type.__setattr__(self, k, v)

    def macro(self):
        for i in self.models:
            i.run()

    def ftosolve(self, x):
        for i, pointer, weight in self._eff_vars:
            self.model[pointer] = x[i] / weight
        self.macro()
        y = [eval(code, model) * weight for code, model, weight in self._eff_cons]
        y = array(y)
        return y

    def ftomin(self, x):
        y = self.ftosolve()
        return sum(y**2)

    def getstart(self):
        x = [v.model[v.pointer] * v.weight for v in self.vars if v.active]
        effv = [(v.pointer, v.weight) for v in self.vars if v.active]
        self._eff_vars = [(i, m, w) for i, (m, w) in enumerate(effv)]
        effc = []
        for v in effv:
            effc.extend(v.getcons())
        effc.extend(self.cons)
        self._eff_cons = [(c.code, c.model, c.weight) for c in effc if c.active]
        self._cons_mask = [c.equality for c in self.cons if c.active]
        return array(x)

    def run(self):
        jacobian(
            self.ftosolve,
            self.getstart(),
            maskstart=self._cons_mask,
            maxsteps=20,
            bisec=4,
            tol=1e-20,
            maxcalls=10000,
            eps=None,
            debug=True,
        )
