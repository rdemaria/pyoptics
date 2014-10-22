from pyviews import *

print 'start'
self=sddsdata('LHCBPM1-Mon_Sep_15_13-51-42_CEST_2008.sdds.sdds.gz','big')
#self=sddsdata('asciisdds.sdds')
#self=sddsdata('binary.sdds.gz','big')
#t=sddsdata('test.sdds')


# x: pad byte (no data);
# c:char;
# b:signed byte;
# B:unsigned byte;
# h:short;
# H:unsigned short;
# i:int;
# I:unsigned int;
# l:long;
# L:unsigned long;
# f:float;
# d:double.
# s:string (array of char);
# p:pascal string (with count byte).
# P:an integer type that is wide enough to hold a pointer.
# q:long long; Q:unsigned long long

#import struct
#def d(f,s):
#  return struct.unpack(f,s[:struct.calcsize(f)])

#fmt='>QQll1003f';struct.unpack(fmt,self.test[:struct.calcsize(fmt)])
