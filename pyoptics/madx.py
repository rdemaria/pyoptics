import cpymad.madx
import re
import numpy as np
from collections import ChainMap


def _to_str(arr, digits, fixed="g"):
    """covert array to string repr"""
    if arr.dtype.kind in "SU":
        return arr
    else:
        fmt = "%%.%d%s" % (digits, fixed)
        out = [fmt % nn for nn in arr]
        return np.array(out)


class Loc:
    def __init__(self, table):
        self.table = table

    def __getitem__(self, key):
        """
        t.loc[1] -> row
        t.loc['a'] -> pattern
        l.loc[10:20]-> range
        l.loc['a':'b'] -> name range
        l.loc['a':'b':'myname'] -> name range with 'myname' column
        l.loc[-2:2:'x'] -> name value range with 'x' column
        l.loc[-2:2:'x',...] -> & combinations
        """
        mask = np.zeros(self.table._nrows, dtype=bool)
        if isinstance(key, int):
            mask[key] = True
        elif isinstance(key, str):
            col = self.table.name
            r = re.compile(key)
            mask[:] = [r.match(nn) for nn in col]
        elif isinstance(key, slice):
            ia = key.start
            ib = key.stop
            ic = key.step
            col = None
            if isinstance(ic, str):
                col = self.table[ic]
            if isinstance(ia, str) or isinstance(ib, str):
                if col is None:
                    col = self.table.name
                if ia is None:
                    ia = 0
                else:
                    r1 = re.compile(ia)
                    ia = np.where([r1.match(nn) for nn in col])[0][0]
                if ib is None:
                    ib = 0
                else:
                    r2 = re.compile(ib)
                    ib = np.where([r2.match(nn) for nn in col])[0][0]
            elif col is None:
                mask[ia:ib:ic] = True
            else:
                if ia is None and ib is None:
                    mask[:] = True
                elif ia is not None and ib is None:
                    mask = col <= ib
                elif ib is not None and ia is None:
                    mask = col >= ia
                else:
                    mask = (col >= ia) & (col <= ib)
        elif isinstance(key, tuple):
            mask = self[key[0]]
            if len(key) > 1:
                mask &= self[key[1:]]
        return mask


class Table(cpymad.madx.Table):
    @property
    def loc(self):
        return Loc(self)

    @property
    def _nrows(self):
        return len(self[self.col_names()[0]])

    def __getitem__(self, column):
        """Get the column data."""
        if isinstance(column, int):
            return self.row(column)
        try:
            return self._cache[column.lower()]
        except KeyError:
            return self.reload(column)

    def eval(self, expr):
        lcl = ChainMap(self, np.__dict__)
        return eval(expr, {}, lcl)

    __call__ = eval

    def mask_name(self, regexp, column="name"):
        """Return a mask for row names matching the given regular expression."""
        regexp = re.compile(regexp)
        idx = np.array([regexp.match(nn) for nn in self[column]], dtype=bool)
        return idx

    def mask_name_range(self, a, b, column="name"):
        r1 = re.compile(a)
        r2 = re.compile(b)
        mask = np.zeros(len(self.name), dtype=bool)
        ia = np.where([r1.match(nn) for nn in self.name])[0][0]
        ib = np.where([r2.match(nn) for nn in self.name])[0][0]
        mask[ia : ib + 1] = True
        return mask

    def mask_value_range(self, a, b, column="s"):
        return np.array([a <= nn <= b for nn in self.eval(column)], dtype=bool)

    def show(
        self,
        rows=None,
        cols=None,
        maxrows=20,
        maxwidth=80,
        output=None,
        digits=6,
        fixed="g",
    ):
        """Pretty print a twiss table

        rows:  `None` all columns,
               a regexp to matching t["name"]-
               a boolean index
        cols:  a list of spaced columns names or expressions
        digits: the number of significant digit to show
        maxrows: maximum number of rows show, None for all
        maxwidth: maximum number row length
        output: None-> stdout, str -> a string, "filename" a file, fh an open file

        Examples:

        tw.show()
        tw.show('mb', 'betx dx/sqrt(betx)')
        tw.show(tw.loc[mb, 10:20:'s'], 'betx dx/sqrt(betx)')
        tw.show(cols='betx', maxwidth=150)
        tw.show(cols='betx', maxrows=None)
        tw.show(cols='betx', digits=12,fixed='e')
        tw.show(output=dict)
        tw.show(output=str)
        tw.show(output=pandas.DataFrame)
        tw.show(output='outfile.txt')
        """
        if rows is None:
            idx = slice(None)
        elif isinstance(rows, str):
            if "name" in self:
                regex = re.compile(rows)
                idx = np.where([regex.match(nn) for nn in self.name])[0]
            else:
                raise (f"Table does not have 'name' column to search")
        elif isinstance(rows, tuple):
            idx = np.where(self.crange(rows[0], rows[1]))[0]
        else:
            rows = np.array(rows)
            if rows.dtype.kind == "b":
                idx = np.where(rows)[0]
            else:
                idx = rows

        cut = -1
        if output is not dict or not hasattr(output, "from_dict"):
            if rows is None and output is not dict and len(self) > maxrows:
                cut = maxrows // 2
                idx = np.r_[np.arange(cut), np.arange(len(self) - cut, len(self))]
            elif hasattr(idx, "__len__") and len(idx) > maxrows:
                cut = maxrows // 2
                idx = np.r_[idx[:cut], idx[-cut:]]

        if cols is None:
            cols = self.col_names()
        elif hasattr(cols, "split"):
            cols = cols.split()

        if "name" not in cols and "name" in self.col_names():
            cols.insert(0, "name")

        if maxwidth is None:
            maxwidth = 1e30

        data = []
        width = 0
        fmt = []
        header = ""
        if output == dict:
            return {cc: self.eval(cc)[idx] for cc in cols}
        if hasattr(output, "from_dict"):
            dct = {cc: self.eval(cc)[idx] for cc in cols}
            dct = output.from_dict(dct)
            if "name" in self and hasattr(dct,'set_index'):
                dct.set_index(self.name)
            return dct

        for cc in cols:
            coldata = self.eval(cc)[idx]
            coltype = coldata.dtype.kind
            col = _to_str(coldata, digits, fixed)
            colwidth = int(col.dtype.str[2:])
            if len(cc) > colwidth:
                colwidth = len(cc)
            colwidth += 1
            width += colwidth + 1
            if width < maxwidth:
                if coltype in "SU":
                    fmt.append("%%-%ds " % (colwidth - 1))
                else:
                    fmt.append("%%%ds" % colwidth)
                header += fmt[-1] % cc
                data.append(col)

        result = [header]
        for ii in range(len(col)):
            row = "".join([ff % col[ii] for ff, col in zip(fmt, data)])
            result.append(row)
            if ii == cut:
                result.append("...")
        result = "\n".join(result)
        if output is None:
            print(result)
        elif output is str:
            return result
        elif hasattr(output, "write"):
            output.write(result)
        else:
            output = pathlib.Path(output)
            with open(output, "w") as fh:
                fh.write(result)


class TableMap(cpymad.madx.TableMap):
    def __getitem__(self, name):
        try:
            return Table(name, self._libmadx)
        except ValueError:
            raise KeyError("Table not found {!r}".format(name)) from None


class Madx(cpymad.madx.Madx):
    def __init__(
        self,
        libmadx=None,
        command_log=None,
        stdout=None,
        history=None,
        prompt=None,
        **Popen_args,
    ):
        super().__init__(libmadx, command_log, stdout, history, prompt, **Popen_args)

        self.table = TableMap(self._libmadx)
