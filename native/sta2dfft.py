#/usr/bin/env python3
#
# JM: 12 Apr 2018
#
# the sta2dfft.f90 adapted for python
# contains 2d spectral commands which uses stafft

from stafft import *

# This module performs FFTs in two directions on two dimensional arrays using 
# the stafft library module to actually compute the FFTs. If FFTs in one 
# direction only are required use the stafft module directly. The module can 
# compute any combination of sine, cosine and full FFTs in each direction. 
# Along with the usual forwards (physical -> Fourier space) and reverse 
# (Fourier space -> physical) routines there are also routines for computing 
# the first derivatives in either direction.
#
# The convention is that for each direction the array is dimensioned 1:nx or 
# 1:ny for either the sine or full transforms. While the cosine transforms 
# require the additional endpoint so 0:nx or 0:ny.
#
# The routines contained in this module are:
#
# init2dfft(nx,ny,lx,ly,xfactors,yfactors,xtrig,ytrig,kx,ky)
#          This routine initialises all the arrays needed for further 
#          transforms. The integers nx and ny are the array dimensions. Then 
#          lx and ly are the domain lengths - these are needed for the correct 
#          scaling when computing derivatives. The arrays xfactors, yfactors, 
#          xtrig and ytrig are needed to perform the various FFTs by the stafft 
#          module (see there for further details. kx and ky are arrays to hold 
#          the wavenumbers associated with each mode in the domain, and are 
#          used in computing derivatives. 
#
#          **If it is known at initialisation that no derivatives are required 
#            it is possible just to pass 1.d0 for each of lx and ly, along with
#            dummy arrays for kx and ky since these are only needed for 
#            computing the derviatives.**

from numpy import pi, arange

#=====================================================================
def init2dfft(nx, ny, lx, ly):
  """
  This subroutine performs the initialisation work for all subsequent
  transform and derivative routines. 
  It calls the initfft() routine from the supproting 1d FFT module for 
  transforms in both x and y directions.
  The routine then defines the two wavenumber arrays, one in each direction.
  
  Input:
    n       = 
  
  Returns:
    factors = 
    
  """
  xfactors, xtrig = initfft(nx)
  yfactors, ytrig = initfft(ny)
  
  if (lx != 0.0):
  # Define x wavenumbers:
    sc = pi / lx
    kx = sc * arange(1, nx + 1)
  else:
  # Catastrophic end to run if wave number definition fails:
    print('**********************************************')
    print(' Wavenumber array definition not possible.')
    print(' Domain length in x equal to zero not allowed.')
    print(' STOPPING...')
    print('**********************************************')

  if (ly != 0.0):
  # Define y wavenumbers:
    sc = pi / ly
    ky = sc * arange(1, ny + 1)
  else:
  # Catastrophic end to run if wave number definition fails:
    print('**********************************************')
    print(' Wavenumber array definition not possible.')
    print(' Domain length in y equal to zero not allowed.')
    print(' STOPPING...')
    print('**********************************************')
    
  return (xfactors, yfactors, xtrig, ytrig, kx, ky)


