class view(object):
    __slots__ = ["data", "index"]

    def __init__(self, data, index=None):
        """ """
        self.data = data
        self.index = index

    def _update(self, v, *args):
        self.data[self.index] = v

    def _env(self, *argl, **argd):
        pass

    def _call(self, *args):
        return self.data[self.index]

    def __repr__(self):
        return "view(...,%s) -> %s" % (repr(self.index), self._call())

    __str__ = __repr__


if __name__ == "__main__":
    import doctest

    doctest.testmod()
