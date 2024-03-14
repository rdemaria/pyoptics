from expr import expr
from view import view


class frame(object):
    """Container class providing customizable update method and deferred
    evaluation. Custom behavior in case of setattr and getattr is called for
    objects providing __call__, _update and _env interfaces.
    - _update is called if an attribute is already in the class
    - _env is called to comunicate to the object the mother class
    - __call__ is called to evaluate the attribute
    Example of implemtation are in view and expr.
    """

    def __init__(self, **args):
        """
        >>> data={1:5}
        >>> a=frame(one=1,two=expr('one+1'),three=view(data,1))
        >>> a.one=1000
        >>> a.two
        1001
        >>> a.three=4
        >>> data[1]
        4
        >>> del a.three
        >>> a.three=expr('two')
        >>> data[1]
        4
        """
        sattr = object.__setattr__
        sattr(self, "_onset", [])
        for k, v in args.items():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        old = self.__dict__.get(k)
        if hasattr(old, "_update"):
            old._update(v)
        else:
            if hasattr(v, "_call"):
                v._env(self)
            object.__setattr__(self, k, v)
        _onset = object.__getattribute__(self, "_onset")
        for i in _onset:
            if hasattr(i, "_call"):
                i._call()
            elif callable(i):
                i()

    def __getattribute__(self, k):
        v = object.__getattribute__(self, k)
        if hasattr(v, "_call"):
            return v._call(self)
        else:
            return v

    def on_set(self, v, env=None):
        if not env:
            env = self
        if hasattr(v, "_call"):
            v._env(self)
            self._onset.append(v)
        for i in self.__dict__.values():
            if isinstance(i, frame):
                i.on_set(v, env)

    def __getitem__(self, k):
        try:
            return getattr(self, k)
        except AttributeError:
            raise KeyError

    def _update(self, v):
        for i in self.__dict__.values():
            if hasattr(i, "_update"):
                i._update(v)

    def _eval(self, expr):
        return eval(expr, locals(), self)

    def __repr__(self):
        out = []
        for k in sorted(self.__dict__):
            if not k.startswith("_"):
                l = []
                l.append(k)
                l.extend(repr(self.__dict__[k]).split("\n"))
                fr = [": ".join(l[:2])]
                fr.extend(["  " + i for i in l[2:]])
                out.append("\n".join(fr))
        return "\n".join(out)

    def _overattr(self, k, v):
        if hasattr(self, k):
            delattr(self, k)
        setattr(self, k, v)

    __str__ = __repr__


if __name__ == "__main__":
    import doctest

    doctest.testmod()
