#/usr/bin/env python3
#
# JM: 11 Apr 2018
#
# just a python script to test the running of stuff

from stafft import *

factors, trig = initfft(128)

print(factors)
for i in range(len(trig)):
  print("%.10f" % trig[i])

