#!/usr/bin/env python3
#
# JM, 29 Oct 2018
#
# some subfunctions to do analysis of pydra data

from pydra_misc import * # numpy is loaded here as np
from casl import parameters, spectral, constants

#-------------------------------------------------------------------------------
# generate the EKE in layers and modes

def calc_eke(data_dir, parameters, constants, kt):

  """
  Subfunction to generate EKE in layers and modes (using zonal mean)
  
  Input:
    data_dir    data directory
    parameters  parameter module from load
    constants   constants module from load
    kt          time stamp
  
  Output:
    eke_L1L2    2d field of EKE in layers
    eke_btbc    2d field of EKE in modes
  """
  
  t_now, qq = read_qq(data_dir, parameters.nx, parameters.ny, kt)
  # swap axis to have the indexing consistent with the Fortran code main_invert
  qq = np.swapaxes(qq, 0, 1)

  # no topography
  fhb = np.zeros((parameters.ny + 1, parameters.nx))

  # invert for the velocity fields
  uu, vv, _ = spectral.main_invert(qq, fhb)
  
  uu_btbc = layers_to_modes(uu, constants)
  vv_btbc = layers_to_modes(vv, constants)
  
  # calculate the EKE
  eke_L1L2 = np.zeros(uu.shape)
  eke_btbc = np.zeros(uu_btbc.shape)
  eke_L1L2[:, :, 0] = zonal_eke(uu[:, :, 0], vv[:, :, 0])
  eke_L1L2[:, :, 1] = zonal_eke(uu[:, :, 1], vv[:, :, 1])
  eke_btbc[:, :, 0] = zonal_eke(uu_btbc[:, :, 0], vv_btbc[:, :, 0])
  eke_btbc[:, :, 1] = zonal_eke(uu_btbc[:, :, 1], vv_btbc[:, :, 1])

  return (eke_L1L2, eke_btbc)
  
#-------------------------------------------------------------------------------
# generate the M and N

def calc_geom_momentum(data_dir, parameters, constants, kt):

  """
  Subfunction to generate M = <v'v' - u'u'> / 2 and N = <u'v'> 
  (using zonal mean)
  
  Input:
    data_dir    data directory
    parameters  parameter module from load
    constants   constants module from load
    kt          time stamp
  
  Output:
    M/N_L1L2    2d field of M/N in layers
    M/N_btbc    2d field of M/N in modes
  """
  
  t_now, qq = read_qq(data_dir, parameters.nx, parameters.ny, kt)
  # swap axis to have the indexing consistent with the Fortran code main_invert
  qq = np.swapaxes(qq, 0, 1)

  # no topography
  fhb = np.zeros((parameters.ny + 1, parameters.nx))

  # invert for the velocity fields
  uu, vv, _ = spectral.main_invert(qq, fhb)
  
  uu_btbc = layers_to_modes(uu, constants)
  vv_btbc = layers_to_modes(vv, constants)
  
  # calculate the quantities
  M_L1L2 = np.zeros(uu.shape)
  N_L1L2 = np.zeros(uu.shape)
  M_btbc = np.zeros(uu_btbc.shape)
  N_btbc = np.zeros(uu_btbc.shape)
  
  zonal_demean(vv[:, :, 0])
  
  M_L1L2[:, :, 0] = (zonal_corr(vv[:, :, 0], vv[:, :, 0])
                   - zonal_corr(uu[:, :, 0], uu[:, :, 0]))
  M_L1L2[:, :, 1] = (zonal_corr(vv[:, :, 1], vv[:, :, 1])
                   - zonal_corr(uu[:, :, 1], uu[:, :, 1]))
  N_L1L2[:, :, 0] =  zonal_corr(vv[:, :, 0], uu[:, :, 0])
  N_L1L2[:, :, 1] =  zonal_corr(vv[:, :, 1], uu[:, :, 1])
  
  M_btbc[:, :, 0] = (zonal_corr(vv_btbc[:, :, 0], vv_btbc[:, :, 0])
                   - zonal_corr(uu_btbc[:, :, 0], uu_btbc[:, :, 0]))
  M_btbc[:, :, 1] = (zonal_corr(vv_btbc[:, :, 1], vv_btbc[:, :, 1])
                   - zonal_corr(uu_btbc[:, :, 1], uu_btbc[:, :, 1]))
  N_btbc[:, :, 0] =  zonal_corr(vv_btbc[:, :, 0], uu_btbc[:, :, 0])
  N_btbc[:, :, 1] =  zonal_corr(vv_btbc[:, :, 1], uu_btbc[:, :, 1])

  return (M_L1L2, N_L1L2, M_btbc, N_btbc)







