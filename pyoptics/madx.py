import cpymad.madx
import re
import numpy as np
from collections import ChainMap

from .tablemixin import TableMixIn

class Table(cpymad.madx.Table,TableMixIn):
    pass


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
