#/usr/bin/env python3
#
# JM: 12 Apr 2018
#
# the secptral adapted for python
# contains inversion commands etc.

from constants import *
# from generic import *
from sta2dfft import *

from numpy import arange, zeros, pi, exp, sqrt, around, log10

# define some parameters
nwx = int(nx / 2)
nwxm1, mwxp1 = nwx - 1, nwx + 1

#========================================================================!
# From main code: call init_invert                   to initialise       !
# then            call main_invert(qq,fhb,uu,vv,pp)  to perform inversion!
#========================================================================!

def init_spectral():
  # set up FFTs:
  xfactors, yfactors, xtrig, ytrig, hrkx, rky = init2dfft(nx, ny, ellx, elly)
  
  
  # Fractional y grid values:
  fac = 1.0 / ny
  yh1 = fac * arange(ny + 1)
  yh0 = 1.0 - yh1
  
  # Define y & beta*y:
  yg = ymin + gly * arange(ny + 1)
  bety = beta * yg

  # Define x wavenumbers:
  rkx = zeros(nx) # define the zero wavenumber here
  for kx in range(1, nwx):
    kxc = nx - kx
    rkx[kx ] = hrkx[2 * kx - 1]
    rkx[kxc] = hrkx[2 * kx - 1]
  rkx[nwx] = hrkx[nx - 1]

  scx = 2.0 * pi / ellx
  rkxmax = scx * nwx
  frkx = zeros(nx) # define the zero here
  wratx = rkx / rkxmax
  frkx[1:nx] = rkx[1:nx] * exp(-36.0 * wratx[1:nx] ** 36.0)

  # Define y wavenumbers:
  scy = pi / elly
  rkymax = scy * ny
  frky = zeros(ny)
  wraty = rky / rkymax
  frky[1:ny] = rky[1:ny] * exp(-36.0 * wraty[1:ny] ** 36.0)
  
  #----------------------------------------------------------------------
  # Initialise arrays for computing the spectrum of any field:
  delk  = sqrt(scx ** 2 + scy ** 2)
  delki = 1.0 / delk
  kmax  = around( sqrt(rkxmax ** 2 + rkymax ** 2) * delki )
  spmf  = zeros(max(nx + 1, ny + 1))
  kmag  = zeros((nx, ny))
  alk   = zeros(max(nx, ny))
  for kx in range(nx):
    k = int( around(rkx[kx] * delki) )
    kmag[kx, 0] = k
    spmf[k] += 1
    
  for ky in range(ny):
    for kx in range(nx):
      k = int(  around( sqrt(rkx[kx] ** 2 + rky[ky] ** 2) * delki )  )
      kmag[kx, ky] = k
      spmf[k] += 1

#  for i in range(len(spmf)):
#    print("%.10f" % spmf[i])
#    print("%13.2f" % spmf[i])
      
  # Compute spectrum multiplication factor (spmf) to account for unevenly
  # sampled shells and normalise spectra by 8/(nx*ny) so that the sum
  # of the spectrum is equal to the L2 norm of the original field:
  snorm = 4.0 * pi / (nx * ny)
  spmf[0] = 0.0
  for k in range(1, int(kmax + 1)):
    spmf[k] = snorm * k / spmf[k]
    alk[k]  = log10(delk * k)
      
  # Only output shells which are fully occupied (k <= kmaxred):
  kmaxred = around(  sqrt( (rkxmax ** 2 + rkymax ** 2) / 2.0 ) *delki  )

#!----------------------------------------------------------------------
# !Define inverse spectral inversion operators:
#if (barot) then
#   !The first mode is barotropic (kd1 = 0):
#  qgop1(0,0)=zero
#else
#  qgop1(0,0)=-one/(kd1sq+small)
#endif
#qgop2(0,0)=-one/kd2sq
#laplace(0,0)=zero
#do kx=1,nxm1
#  rksq=rkx(kx)**2
#  qgop1(kx,0)=-one/(rksq+kd1sq)
#  qgop2(kx,0)=-one/(rksq+kd2sq)
#  laplace(kx,0)=-rksq
#enddo
#do ky=1,ny
#  rksq=rky(ky)**2
#  qgop1(0,ky)=-one/(rksq+kd1sq)
#  qgop2(0,ky)=-one/(rksq+kd2sq)
#  laplace(0,ky)=-rksq
#enddo
#do ky=1,ny
#  do kx=1,nxm1
#    rksq=rkx(kx)**2+rky(ky)**2
#    qgop1(kx,ky)=-one/(rksq+kd1sq)
#    qgop2(kx,ky)=-one/(rksq+kd2sq)
#    laplace(kx,ky)=-rksq
#  enddo
#enddo

#!----------------------------------------------------------------------
# !Hyperbolic functions used to correct boundary conditions in inversion:
# !First mode:
#do kx=1,nxm1
#  fac=sqrt(rkx(kx)**2+kd1sq)*elly
#  div=one/(one-exp(-two*fac))
#  do iy=1,nym1
#    argm=fac*(one-yh1(iy))
#    argp=fac*(one+yh1(iy))
#    decy1(iy,kx)=(exp(-argm)-exp(-argp))*div
#  enddo
#enddo

# !Second mode (kd2 > kd1 is assumed):
#do kx=1,nxm1
#  fac=sqrt(rkx(kx)**2+kd2sq)*elly
#  div=one/(one-exp(-two*fac))
#  do iy=1,nym1
#    argm=fac*(one-yh1(iy))
#    argp=fac*(one+yh1(iy))
#    decy2(iy,kx)=(exp(-argm)-exp(-argp))*div
#  enddo
#enddo


