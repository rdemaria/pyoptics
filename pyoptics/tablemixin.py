import re
import numpy as np
from collections import ChainMap
import pathlib


def _to_str(arr, digits, fixed="g"):
    """covert array to string repr"""
    if arr.dtype.kind in "SU":
        return arr
    elif arr.dtype.kind == "O":
        return arr.astype("U")
    elif arr.dtype.kind in "iu":
        return np.char.mod("%d", arr)
    else:
        fmt = "%%.%d%s" % (digits, fixed)
        return np.char.mod(fmt, arr)


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
        elif hasattr(key,'dtype'):
            mask=key
        elif isinstance(key, str):
            col = self.table._index
            r = re.compile(key, re.IGNORECASE)
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
                    col = self.table._index
                if ia is not None:
                    r1 = re.compile(ia, re.IGNORECASE)
                    ia = np.where([r1.match(nn) for nn in col])[0][0]
                if ib is not None:
                    r2 = re.compile(ib, re.IGNORECASE)
                    ib = np.where([r2.match(nn) for nn in col])[0][0] + 1
                mask[ia:ib] = True
            elif col is None:
                mask[ia:ib:ic] = True
            else:
                if ia is None and ib is None:
                    mask |= True
                elif ia is not None and ib is None:
                    mask |= col <= ib
                elif ib is not None and ia is None:
                    mask |= col >= ia
                else:
                    mask |= (col >= ia) & (col <= ib)
        elif isinstance(key, tuple):
            mask = self[key[0]]
            if len(key) > 1:
                mask &= self[key[1:]]
        return mask


class TableMixIn:
    """
    Assumptions:
    - __getitem__[cname]: return a column
    - __contains__[cname]: return if column exists
    - col_names(): return list of column names
    - row: return a dictionary at row
    - _cached: dictionary cointaining data
    - reload(name): return array
    """

    @property
    def loc(self):
        return Loc(self)

    @property
    def _nrows(self):
        return len(self[self.col_names()[0]])

    @property
    def _index(self):
        if not hasattr(self, "_index_name"):
            if "name" in self:
                return self["name"]
            else:
                return self[self.col_names[0]]
        else:
            return self[self._index_name]

    def set_index(self, name):
        self._index_name = name

    def __getitem__(self, column):
        """Get the column data."""
        if isinstance(column, int):
            return self.row(column)
        try:
            return self._cache[column.lower()]
        except KeyError:
            if hasattr(self, column):
                return getattr(self, column)
            else:
                return self.reload(column)

    def eval(self, expr):
        lcl = ChainMap(self, np.__dict__)
        return eval(expr, {}, lcl)

    __call__ = eval

    def _idx_from_regex(self, regex, index=None):
        return np.where(self.mask(regex, index))[0]

    def mask(self, regex, index=None):
        if index is None:
            col = self._index
        else:
            col = self[index]
        regex = re.compile(regex, re.IGNORECASE)
        return np.fromiter((regex.match(nn) for nn in col), dtype="bool")

    __floordiv__ = mask

    def __len__(self):
        return self._nrows

    def show(
        self,
        rows=None,
        cols=None,
        maxrows=20,
        maxwidth=None,
        output=None,
        digits=6,
        fixed="g",
    ):
        """Pretty print a twiss table

        rows:  `None` all columns,
               a regexp to matching t["name"]
               a boolean index
        cols:  a list of spaced columns names or expressions
        digits: the number of significant digit to show
        maxrows: maximum number of rows show, None for all
        maxwidth: maximum number row length
        output: None-> stdout, str -> a string, "filename" a file, fh an open file
        fixed: 'g' or 'f'

        Examples:

        tw.show()
        tw.show('mb', 'betx dx/sqrt(betx)')
        tw.show(tw.loc['mb', 10:20:'s'], 'betx dx/sqrt(betx)')
        tw.show(cols='betx', maxwidth=150)
        tw.show(cols='betx', maxrows=None)
        tw.show(cols='betx', digits=12,fixed='e')
        tw.show(output=dict)
        tw.show(output=str)
        tw.show(output=pandas.DataFrame)
        tw.show(output='outfile.txt')
        """

        tlen = self._nrows

        index_name = self._index_name if hasattr(self, "_index_name") else "name"

        if cols is None and maxwidth is None:
            maxwidth = 80

        if rows is None:
            idx = slice(None)
        elif isinstance(rows, str):
            idx = self._idx_from_regex(rows)
        elif isinstance(rows, tuple):
            idx = np.where(self.loc[rows[0] : rows[1]])[0]
        else:
            rows = np.array(rows)
            if rows.dtype.kind == "b":
                idx = np.where(rows)[0]
            else:
                idx = rows

        cut = -1
        if (
            maxrows is not None
            and output is not dict
            and not hasattr(output, "from_dict")
        ):
            if rows is None and output is not dict and tlen > maxrows:
                cut = maxrows // 2
                idx = np.r_[np.arange(cut), np.arange(tlen - cut, tlen)]
            elif hasattr(idx, "__len__") and len(idx) > maxrows:
                cut = maxrows // 2
                idx = np.r_[idx[:cut], idx[-cut:]]

        if cols is None:
            cols = self.col_names()
        elif hasattr(cols, "split"):
            cols = cols.split()

        if index_name not in cols and index_name in self.col_names():
            cols.insert(0, "name")

        if maxwidth is None:
            maxwidth = 1e30

        data = []
        width = 0
        fmt = []
        header = []
        if output == dict:
            return {cc: self.eval(cc)[idx] for cc in cols}
        if hasattr(output, "from_dict"):
            dct = {cc: self.eval(cc)[idx] for cc in cols}
            dct = output.from_dict(dct)
            if index_name in self and hasattr(dct, "set_index"):
                dct.set_index(index_name)  # dataframe
            return dct

        for cc in cols:
            coldata = self.eval(cc)[idx]
            coltype = coldata.dtype.kind
            col = _to_str(coldata, digits, fixed)
            colwidth = int(col.dtype.str[2:])
            if len(cc) > colwidth:
                colwidth = len(cc)
            width += colwidth + 1
            if width < maxwidth:
                if coltype in "SU":
                    fmt.append("%%-%ds" % (colwidth))
                else:
                    fmt.append("%%%ds" % colwidth)
                header.append(fmt[-1] % cc)
                data.append(col)

        result = [" ".join(header)]
        for ii in range(len(col)):
            row = " ".join([ff % col[ii] for ff, col in zip(fmt, data)])
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
