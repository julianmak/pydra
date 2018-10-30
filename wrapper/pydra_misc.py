#!/usr/bin/env python3
#
# JM, 20 Apr 2018
#
# some misc functions relating to pydra

import os
import numpy as np
from numpy.fft import rfft, irfft

#-------------------------------------------------------------------------------
# read the PV data
def read_qq(data_dir, nx, ny, kt, num_frame = False):

  """
  Subfunction to read the PV data (assumes they are called qq1 and qq2)
  
  Input:
    data_dir    string for where data is located
    nx          number of x points
    ny          number of y points
    kt          index of frame to pull out (set to 0 to see how many frames are found)
  
  Output:
    t_now            actual dump time from header
    qq[x, y, layer]  PV data (NOTE: indexing may change)
  """
  
  N = nx * (ny + 1) + 1
  
  qq = np.zeros((nx, ny + 1, 2))

  qq1_filename = data_dir + "qq1.r4"

  qq1_file = open(qq1_filename, "r")
  file_bytes = os.path.getsize(qq1_filename)
  nframes = int(file_bytes / (4 * N))

  if num_frame:
    print("number of frames found = %i " % nframes)

  raw_array = np.fromfile(qq1_file, dtype = np.float32)
  time_eles = [i * N for i in range(nframes)]

  # pull out the header time
  if kt == nframes - 1:
    index_start, index_end = time_eles[kt] + 1, len(raw_array)
  else:
    index_start, index_end = time_eles[kt] + 1, time_eles[kt + 1]
  t_now = raw_array[time_eles[kt]]

  # pull out the data relevant to time and reshape
  qq[:, :, 0] = raw_array[index_start:index_end].reshape((nx, ny + 1))
  
  qq1_file.close()
  
  qq2_filename = data_dir + "qq2.r4"
  qq2_file = open(qq2_filename, "r")
  raw_array = np.fromfile(qq2_file, dtype = np.float32)
  qq[:, :, 1] = raw_array[index_start:index_end].reshape((nx, ny + 1))
  
  qq2_file.close()
  
  return (t_now, qq)
  
#-------------------------------------------------------------------------------
# transform between layer fields to modal fields

def layers_to_modes(field_L1L2, constants):

  """
  Subfunction to transform layer fields to modal fields
  
  Input:
    field_L1L2 2d field in layers
  
  Output:
    field_btbc 2d field in modes
  """
  
  field_btbc = np.zeros(field_L1L2.shape)
  
  field_btbc[:, :, 0] = (constants.vec11 * field_L1L2[:, :, 0] 
                       + constants.vec12 * field_L1L2[:, :, 1])
  field_btbc[:, :, 1] = (constants.vec21 * field_L1L2[:, :, 0] 
                       + constants.vec22 * field_L1L2[:, :, 1])

  return field_btbc
  
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
  
  outfield = np.mean(infield, axis = 1)
  
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
  outfield = np.zeros(infield.shape)
  
  if infield.shape[1] < 10:
    print("WARNING: x-axis should be longer than this, check you are throwing in arrays of the right shape!")
  
  for i in range(infield.shape[1]):
    outfield[:, i] = infield[:, i] - mean_1d # probably a faster way of doing this (broadcast?)
  
  return outfield

#-------------------------------------------------------------------------------
def zonal_corr(f1_in, f2_in):
  """
  Subfunction to generate the correlation function of two input fields,
  as defined by the zonal average
  
  e.g. u_in and v_in would give Reynolds stress
       u_in and u_in would give part of M in geometric decomposition
  
  can do (u_in, u_in) and (v_in, v_in) to give eke
  
  Input:
    f1_in, f2_in 2d input field
  
  Output:
    f1f2_out     2d output demean-ed field
  """
  
  # do a de-mean of both fields even if some of them have no mean 
  # (e.g. v = d(psi)/dx)
  f1_fluc = zonal_demean(f1_in)
  f2_fluc = zonal_demean(f2_in)
  
  f1f2_out = f1_fluc * f2_fluc
  
  return f1f2_out
  
#-------------------------------------------------------------------------------
def zonal_eke(u_in, v_in):
  """
  Subfunction to generate the eke field as defined by the zonal average
  
  Input:
    u_in, v_in  2d input velocity field
  
  Output:
    eke         2d output eke field
  """
  
  eke = (zonal_corr(u_in, u_in) + zonal_corr(v_in, v_in)) / 2.0
  
  return eke
  
#-------------------------------------------------------------------------------
def zonal_eke_int(u1_in, v1_in, u2_in, v2_in, parameters):
  """
  Subfunction to compute the domain-averaged eke value, where the
  eke is defined by the zonal average (give it the layer velocity though)
  
  Input:
    u_in, v_in  2d layer input velocity field
    parameters  collection of system parameters
  
  Output:
    eke_domavg  domain integrated value
    (TODO: eke_max?)
  """
  # work out the eke of layers then weight it with layer depth (which is 1)
  eke1 = zonal_eke(u1_in, v1_in) * parameters.h1
  eke2 = zonal_eke(u2_in, v2_in) * (1 - parameters.h1)
  
  # add them together and take a mean (divide by depth = 1 not written in)
  eke_domavg = np.mean(np.mean(eke1 + eke2, axis = 0), axis = 0)
  
  return eke_domavg
  
#-------------------------------------------------------------------------------
def zonal_ens(q_in):
  """
  Subfunction to generate the enstrophy field as defined by the zonal average
  
  Input:
    q_in        2d input PV field
  
  Output:
    ens         2d output enstrophy field
  """
  
  ens = zonal_corr(q_in, q_in) / 2.0
  
  return ens

#-------------------------------------------------------------------------------
def zonal_ens_int(q1_in, q2_in, parameters):
  """
  Subfunction to compute the domain-averaged ens value, where the
  ens is defined by the zonal average (give it the layer PV though)
  
  Input:
    q_in        2d layer input PV field
    parameters  collection of system parameters
  
  Output:
    ens_domavg  domain integrated value
    (TODO: ens_max?)
  """
  # work out the eke of layers then weight it with layer depth (which is 1)
  ens1 = zonal_ens(q1_in) * parameters.h1
  ens2 = zonal_ens(q2_in) * (1 - parameters.h1)
  
  # add them together and take a mean (divide by depth = 1 not written in)
  ens_domavg = np.mean(np.mean(ens1 + ens2, axis = 0), axis = 0)
  
  return ens_domavg
  
#-------------------------------------------------------------------------------
def zonal_mode_extract(infield, mode_keep, low_pass = False):
  """
  Subfunction to extract or swipe out zonal modes (mode_keep) of (y, x) data.
  Assumes here that the data is periodic in axis = 1 (in the x-direction) with
  the end point missing
  
  If mode_keep = 0 then this is just the zonal averaged field
  
  Input:
    in_field    2d layer input field
    mode_keep   the zonal mode of the data to be extracted from
    
  Opt input:
    low_pass    get rid of all modes from mode_keep + 1 onwards
  
  Output:
    outfield    zonal mode of the data
  """
  
  outfield_h = rfft(infield, axis = 1)
  outfield_h[:, mode_keep+1::] = 0
  if not low_pass:
    outfield_h[:, 0:mode_keep] = 0
  
  return irfft(outfield_h, axis = 1)
