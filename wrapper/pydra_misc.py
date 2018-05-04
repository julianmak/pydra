#!/usr/bin/env python3
#
# JM, 20 Apr 2018
#
# some misc functions relating to pydra

import os
from numpy import fromfile, float32, delete, zeros

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
  
  qq = zeros((nx, ny + 1, 2))

  qq1_filename = data_dir + "qq1.r4"

  qq1_file = open(qq1_filename, "r")
  file_bytes = os.path.getsize(qq1_filename)
  nframes = int(file_bytes / (4 * N))

  if num_frame:
    print("number of frames found = %i " % nframes)

  raw_array = fromfile(qq1_file, dtype = float32)
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
  raw_array = fromfile(qq2_file, dtype = float32)
  qq[:, :, 1] = raw_array[index_start:index_end].reshape((nx, ny + 1))
  
  qq2_file.close()
  
  return (t_now, qq)
  
#-------------------------------------------------------------------------------
# transform between layer fields to modal fields

def layers_to_modes(field_top, field_low, constants):

  """
  Subfunction to transform layer fields to modal fields
  
  Input:
    field_top   2d field of top layer
    field_low   2d field of lower layer
  
  Output:
    field_bt    2d field of barotropic mode
    field_bc    2d field of baroclinic mode
  """
  
  field_bt = constants.vec11 * field_top + constants.vec12 * field_low
  field_bc = constants.vec21 * field_top + constants.vec22 * field_low
  
  return (field_bt, field_bc)
