#!/usr/bin/env python3
#
# JM, 20 Apr 2018
#
# some subfunctions to average things

from numpy import zeros, mean

#-------------------------------------------------------------------------------
def zonal_ave(infield):
  """
  Subfunction to zonally average some input 2d field. Assumes it is uniformly 
  spaced in x (axis = 1)
  
  Input:
    infield     2d input field for data
  
  Output:
    outfield    1d zonally averaged data
  """
  
  outfield = mean(infield, axis = 1)
  
  return outfield
  
#-------------------------------------------------------------------------------
def zonal_demean(infield):
  """
  Subfunction to get rid of the zonal mean of the input field
  
  Input:
    infield     2d input field for data
  
  Output:
    outfield    2d output field for data
  """
  
  mean_1d = zonal_ave(infield)
  field_mean = zeros(infield.shape)
  
  for i in range(infield.shape[1]):
    field_mean[:, i] = mean_1d # probably a faster way of doing this (broadcast?)
  
  outfield = infield - field_mean
  
  return outfield
  
  
