module parameters
! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.
 !Initial time loop from which to begin the simulation:
integer,parameter:: loopinit=0
 !loopinit > 0 may be used for continuing a simulation (see replace below)
 !Logical to indicate replacing existing data with the output of a new
 !simulation or to append to existing data:
logical,parameter:: replace=.false.
 !Number of fluid layers (here always 2):
integer,parameter:: nz=2
! Number of grid boxes in the x & y directions (inversion grid):
integer,parameter:: nx=128,ny=64
 !Number of contours used for representing PV in either layer:
integer,parameter:: ncontq=80
! ncontq : used to compute the PV contour interval from
! dq = (qq_max-qq_min)/ncontq (in each layer)
! Simulation time length etc..
double precision,parameter:: tsim=750.0
double precision,parameter:: tgsave=1.0,tcsave=50.0
! tsim : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment (approximate)
! ***NOTE*** tcsave should always be an integer multiple of tgsave
 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816
!***Physical parameters:***
double precision,parameter:: alpha=1.0,h1=0.50,kdbar=20.0
double precision,parameter:: ellx=6.28318530717958647692
double precision,parameter:: elly=3.14159265358979323846
double precision,parameter:: ymin=-1.57079632679489661923
double precision,parameter:: beta=2.51327412287183459076
double precision,parameter:: u1bot=0.0,u1top=0.0
double precision,parameter:: u2bot=0.0,u2top=0.0
double precision,parameter:: rtherm1=0.0
double precision,parameter:: rtherm2=0.0
double precision,parameter:: rekman=0.0
double precision,parameter:: eirate=0.0
double precision,parameter:: vorvor=1.0
double precision,parameter:: rheton=0.1
double precision,parameter:: thdmax=0.0
double precision,parameter:: sponlen=3.14159265358979323846
integer,parameter:: iseed=2583211,mod1bc=0
logical,parameter:: topogr=.false.
! h1 : fractional thickness of the lower layer
! kdbar : f*sqrt{(H1+H2)/(g*H1*H2*(1-alpha))} where H1 & H2 are the
! layer depths, f is the Coriolis frequency and g is gravity
! ellx : domain width in x (periodic, centred at 0)
! elly : domain width in y
! ymin : minimum value of y (ymax = ymin + elly; see constants.f90)
! beta : planetary vorticity gradient
! u1bot : zonal mean mode 1 (smaller kd) velocity at y = ymin
! u1top : zonal mean mode 1 (smaller kd) velocity at y = ymax
! *** NOTE: these are not used if kd = kd1 = 0 ***
! u2bot : zonal mean mode 2 (larger kd) velocity at y = ymin
! u2top : zonal mean mode 2 (larger kd) velocity at y = ymax
! rtherm1: thermal damping rate (1/tau_1) of lower layer thickness
! rtherm2: thermal damping rate (1/tau_2) of upper layer thickness
! rekman : Ekman damping rate (1/tau_E)
! eirate : enstrophy input rate (via heton pairs of point vortices
! which are converted to gridded vorticity and added to qd)
! vorvor : the mean vorticity magnitude associated with the point
! vortex upon placement on the grid
! rheton : the radius of each of the pair of vortices added as a heton
! iseed : seed for initialising point vortex forcing
! topogr : logical variable to indicate presence of topography
! mod1bc : boundary condition used on mode 1 when kd1 > 0;
! use 1 to implement fixed zonal mean velocities u1bot & u1top;
! use 2 to implement zero average streamfunction and velocity;
! (not used if kd1 = 0)
! thdmax : max rate of thermal damping at the domain edges
! sponlen: sponge length of vorticity damping; the function used
! is defined in casl.f90 (see "sponge")
!----------------------------------------------------------------
end module
