import gzip

import numpy as n


def _myopen(fn):
    try:
        if fn.endswith(".gz"):
            return gzip.open(fn)
        else:
            return open(fn)
    except IOError:
        import io

        return io.StringIO(fn)


yasptypes = {"%s": str, "%d": int, "%f": float}
coltypes = {
    "BEAM": "i",
    "HW-STATUS": "i",
    "KICK": "d",
    "NAME": "S20",
    "PLANE": "S1",
    "POS": "d",
    "RMS": "d",
    "STATUS": "i",
    "STATUS-TAG": "S20",
    "STRENGTH-NAME": "S20",
    "SUM": "d",
    "RT-KICK": "d",
}


def readdata(o, cols, count):
    ds = []
    for i in range(count):
        ds.append(o.readline().split())
    ds = zip(*ds)
    d = {}
    for name, data in zip(cols, ds):
        d[name.lower()] = n.array(data, dtype=coltypes[name])
    return d


class YASPData(object):
    def __init__(self, filename):
        self.filename = filename
        o = _myopen(filename)
        c = ""
        header = {}
        dataset = {}
        l = o.readline()
        while l:
            c = l[0]
            if c == "@":
                #        print l
                l = l.strip().split(None, 4)
                if len(l) == 4:
                    name, name, type, data = l
                elif len(l) == 2:
                    name = l[1]
                    type = "%d"
                    data = 1
                header[name] = yasptypes[type](data)
            elif c == "#":
                name, name, type, data = l.strip().split(None, 4)
                cols = o.readline().strip().split()
                cols.pop(0)
                if name == "MONITOR":
                    dataset["monitor-h"] = readdata(o, cols, header["MONITOR-H-NUM"])
                    dataset["monitor-v"] = readdata(o, cols, header["MONITOR-V-NUM"])
                elif name == "CORRECTOR":
                    dataset["corrector-h"] = readdata(
                        o, cols, header["CORRECTOR-H-NUM"]
                    )
                    dataset["corrector-v"] = readdata(
                        o, cols, header["CORRECTOR-V-NUM"]
                    )
            l = o.readline()
        self.dataset = dataset
        self.header = header
