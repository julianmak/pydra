##/usr/bin/env python3
#
# JM: 11 Apr 2018
#
# constants that are unchanged, depending on the customisable parameters.py
from parameters import *
from numpy import sqrt, pi

# Contains all the non-modifiable parameters as well as all 
# quantities which never change throughout a simulation
# for the suite of casl f90 codes.

# Grid dimensions +/- 1 & 2:
nxp1, nxm1, nxm2 = nx + 1, nx - 1, nx - 2
nyp1, nym1       = ny + 1, ny - 1

# Fine grid used normally in contour -> grid conversion: 
mgf           = 4
nxf, nyf      = mgf * nx, mgf * ny
# mgf:  fine grid/coarse grid ratio (4 is required by subroutine 
#       coarsen in contours.f90)
# Fine grid dimensions +/- 1 & 2:
nxfm1, nxfm2  = nxf - 1, nxf - 2
nyfp1, nyfm1  = nyf + 1, nyf - 1

# Ultra-fine grid used in contouring: 
mgu           = 16
nxu, nyu      = mgu * nx, mgu * ny
# mgu:  ultra-fine grid/coarse grid ratio (16 is the default)
# Ultra-fine grid dimensions +/- 1 & 2:
nxum1, nxum2  = nxu - 1, nxu - 2
nyup1, nyum1  = nyu + 1, nyu - 1

# For reading & writing direct access data:
ngridp        = nx * nyp1
nbytes        = 4 * (ngridp + 1)

# Maximum number of contour levels (used in surgery and congen):
nlevm         = 2000
# nlevm: up to 2*nlevm contour levels are allowed

# Maximum number of contour nodes:
npm           = 625 * nx * ny
# Maximum number of contours:
nm            = npm / 20 + npm / 200
# Maximum number of nodes on any single contour:
nprm          = npm / 10
# Maximum number of nodes in any contour level:
nplm          = npm / 2

# Generic double precision numerical constants: 
# JM: got rid of this

# Domain lengths and inverses:
hlx           = 0.5 * ellx
hlxi          = 1.0 / hlx
xmin, xmax    = -hlx, hlx
ymax          = ymin + elly
ycen          = (ymin + ymax) / 2.0
hly           = (1.0 - 1.0e-12) * 0.5 * elly
ybeg          = ycen - hly

# Fractional thickness of the upper layer (note: h1 + h2 = 1):
h2            = 1.0 - h1
h1h2          = h1 * h2
h1inv, h2inv  = 1.0 / h1, 1.0 / h2

# Constants used to define relative vorticity and interface displacements:
alphac        = 1.0 - alpha
kdbarsq       = kdbar * kdbar
h1kdbarsq     = h1 * kdbarsq
h1h2kdbarsq   = h1h2 * kdbarsq
h1h2ackdbarsq = h1h2 * alphac * kdbarsq
massfac       = 1.0 / (h1 + alpha * h2)

# Define Rossby deformation wavenumbers (for each mode):
gamma1        = 0.50 - sqrt(0.25 - alphac * h1 * h2)
gamma2        = 0.50 + sqrt(0.25 - alphac * h1 * h2)
kd1sq         = gamma1 * kdbarsq
kd1           = sqrt(kd1sq)
kd2sq         = gamma2 * kdbarsq
kd2           = sqrt(kd2sq)

# Define layer -> mode transformation coefficients:
vec11         = h1 / (1.0 - gamma1)            #mode 1, layer 1
vec12         = (h2 - gamma1) / (1.0 - gamma1) #mode 1, layer 2
vec21         = h1 / (h2 - gamma2)             #mode 2, layer 1
vec22         = 1.0                            #mode 2, layer 2
vec2sum       = vec21 + vec22

# Define mode -> layer transformation coefficients:
determinant   = vec11 * vec22 - vec12 * vec21
vect11        = vec22 / determinant       #layer 1, mode 1
vect12        =-vec12 / determinant       #layer 1, mode 2
vect21        =-vec21 / determinant       #layer 2, mode 1
vect22        = vec11 / determinant       #layer 2, mode 2

# Initialise coefficients for decomposing energy & enstrophy into modes:
coeff1        = h1 * (vect11 ** 2) + alpha * h2 * (vect21 ** 2)
coeff2        = h1 * (vect12 ** 2) + alpha * h2 * (vect22 ** 2)

# Maximum time step:
dtmax         = 0.1 * min(tgsave,tcsave)

# Set maximum Rossby wave frequency (used in adapt in evolution.f90):
srwfm         = beta / ( (2.0 * pi / ellx)**2 + kd1sq )

# Basic constants:
domarea       = ellx * elly
aspect        = ellx / elly
glx, glxi     = ellx / nx, nx / ellx
gly, glyi     = elly / ny, ny / elly
garea, dsumi  = glx * gly, 1.0 / (nx * ny)

 #Logical control variables:
sponge        = (thdmax > 0.0)
heating       = ( (rtherm1 > 0.0) or (rtherm2 > 0.0) or sponge)
friction      = (rekman > 0.0)
damping       = (heating or friction)
barot         = (alpha == 1.0)
stoch         = (eirate > 0.0)


