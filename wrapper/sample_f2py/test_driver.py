#!/usr/bin/env python3
#
from numpy import f2py

sourcefile = open("test_mod.f90", "rb") # open as binary
sourcecode = sourcefile.read()
print(sourcecode)

f2py.compile(sourcecode, modulename = "murp", extension = ".f90")

from murp import test_mod

test_mod.print_murp(6.0)

p = test_mod.add_murp(3.0, 4.0)

print(p)

j = test_mod.self_murp(5.0)

print(j)

p, q = test_mod.multiple_murp(6.0, 7.0)

print(p)

print(q)
