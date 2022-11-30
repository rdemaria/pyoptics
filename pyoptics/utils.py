from numpy import log10, pi, sqrt, exp, isreal
import gzip, io
import re


def pythonname(string):
    string = string.replace("[", "")
    string = string.replace("]", "")
    string = string.replace(".", "_")
    string = string.replace("$", "_")
    string = string.lower()
    return string


def eng_string(x, format="%s", si=False):
    import math

    """
    Returns float/int value <x> formatted in a simplified engineering format -
    using an exponent that is a multiple of 3.

    format: printf-style string used to format the value before the exponent.

    si: if true, use SI suffix for exponent, e.g. k instead of e3, n instead of
    e-9 etc.

    E.g. with format='%.2f':
        1.23e-08 => 12.30e-9
             123 => 123.00
          1230.0 => 1.23e3
      -1230000.0 => -1.23e6

    and with si=True:
          1230.0 => 1.23k
      -1230000.0 => -1.23M
    """
    x = float(x)
    sign = ""
    if x == 0:
        return format % x
    if x < 0:
        x = -x
        sign = "-"
    exp = int(math.floor(math.log10(x)))
    exp3 = exp - (exp % 3)
    x3 = x / (10**exp3)

    if si and exp3 >= -24 and exp3 <= 24 and exp3 != 0:
        exp3_text = "yzafpnum kMGTPEZY"[(exp3 - (-24)) // 3]
    elif exp3 == 0:
        exp3_text = ""
    else:
        exp3_text = "e%03d" % exp3

    # print(exp3_text)
    # print( ( '%s'+format+'%s') % ( sign, x3, exp3_text) )
    return ("%s" + format + "%s") % (sign, x3, exp3_text)


def numtostr(n, ns=12, np=3):
    """Convert a number in a string where . has a fixed position
    np minimum precision
    ns string size

    >>> for i in range(-6,10):
    ...   print numtostr( .1234*10**i),numtostr( .1234567*10**i,np=6)
    ...   print numtostr(-.1234*10**i),numtostr(-.1234567*10**i,np=6)
     123.400e-09  123.456700e-09
    -123.400e-09 -123.456700e-09
       1.234e-06    1.234567e-06
      -1.234e-06   -1.234567e-06
      12.340e-06   12.345670e-06
     -12.340e-06  -12.345670e-06
     123.400e-06  123.456700e-06
    -123.400e-06 -123.456700e-06
       1.234e-03    1.234567e-03
      -1.234e-03   -1.234567e-03
      12.340e-03   12.345670e-03
     -12.340e-03  -12.345670e-03
     123.400e-03  123.456700e-03
    -123.400e-03 -123.456700e-03
       1.234        1.234567
      -1.234       -1.234567
      12.340       12.345670
     -12.340      -12.345670
     123.400      123.456700
    -123.400     -123.456700
       1.234e+03    1.234567e+03
      -1.234e+03   -1.234567e+03
      12.340e+03   12.345670e+03
     -12.340e+03  -12.345670e+03
     123.400e+03  123.456700e+03
    -123.400e+03 -123.456700e+03
       1.234e+06    1.234567e+06
      -1.234e+06   -1.234567e+06
      12.340e+06   12.345670e+06
     -12.340e+06  -12.345670e+06
     123.400e+06  123.456700e+06
    -123.400e+06 -123.456700e+06
    """
    return eng_string(n, "%%%d.%df" % (ns - 5, np), si=False)
    n = float(n)
    #  if abs(n)>0:
    #    l=log10(abs(n))
    #    if l<0:
    #      l=int(l)
    #      o=(l-1)//3*3
    #      fmt='%%%d.%dfe%+03d' % (np+5,np,o)
    #      n=n/10**o
    ##    fill space of digits
    ##    elif -3<=l and l< 0: fmt='%%12.%df' % (np+4)
    ##    elif  0<=l and l< 3: fmt='%%12.%df' % (np+4)
    ##    elif -3<=l and l< 0: fmt='%%11.%df ' % (np+3)
    #    elif  0<=l and l< 3: fmt='%%%d.%df    ' % (np+5,np)
    #    elif  3<=l:
    #      l=int(l)
    #      o=(l)//3*3
    #      fmt='%%%d.%dfe%+03d' % (np+5,np,o)
    #      n=n/10**o
    #  else:
    #    fmt='%4.0f.'+' '*(np+4)
    np = min(np, ns - 6)
    if abs(n) == 0:
        return ("%%%d.%%df" % (ns, np)) % n
    else:
        l = log10(abs(n))
        if l >= -4 and l < 5:
            return ("%%%d.%df" % (ns, np)) % n
        else:
            return ("%%%d.%dg" % (ns, np)) % n


def gt(a, b):
    if a < b:
        return (b - a) ** 2
    else:
        return 0.0


def lt(a, b):
    if a > b:
        return (b - a) ** 2
    else:
        return 0.0


def cmp(a, b, c):
    if a < b:
        return (b - a) ** 2
    elif a > c:
        return (a - c) ** 2
    else:
        return 0.0


def eq(a, b):
    return (a - b) ** 2


def rng(a, b, c):
    return (a > b) & (a < c)


def mystr(d, nd, np):
    """truncate a number or a string with a fix number of chars
    >>> for d in [ 0.443,'stre',1.321E-4,-3211233]: print mystr(d,12)
     443.000e-03
    stre
     132.100e-06
      -3.211e+06
    """
    if isreal(d):
        d = numtostr(d, nd, np)
    return ("%%-%d.%ds" % (nd, nd)) % d


def myflatten(lst):
    for elem in lst:
        if type(elem) in (tuple, list):
            for i in myflatten(elem):
                yield i
        else:
            yield elem


import time


def timeit(s, set=False):
    if set:
        timeit.mytime = time.time()
    print("%8.3f %s" % (time.time() - timeit.mytime, s))


timeit.mytime = 0


def myopen(fn):
    try:
        if fn.endswith(".gz"):
            return gzip.open(fn)
        else:
            return open(fn)
    except IOError:
        return io.StringIO(fn)


def no_dots(x):
    return x.group().replace(".", "_")


madname = re.compile(r"([a-z_][a-z_0-9\.]*)")


def pyname(n):
    n = n.lower()
    n = madname.sub(no_dots, n)
    n.replace("^", "**")
    if n == "from":
        n = "From"
    return n


if __name__ == "__main__":
    import doctest

    doctest.testmod()
