module spectral

use parameters
use constants
use generic
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Maximum x wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1

 !Common arrays, constants:
double precision:: yg(0:ny),pbar(0:ny),bety(0:ny)
double precision:: decy1(nym1,0:nxm1),decy2(nym1,0:nxm1)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: pftop1(0:ny),pfbot1(0:ny),uftop1(0:ny),ufbot1(0:ny)
double precision:: pftop2(0:ny),pfbot2(0:ny),uftop2(0:ny),ufbot2(0:ny)
double precision:: qgop1(0:nxm1,0:ny),qgop2(0:nxm1,0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)
double precision:: facp1,facu1

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: frkx(0:nxm1),frky(ny)
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(0:nxm1,0:ny),kmax,kmaxred

!========================================================================!
! From main code: call init_invert                   to initialise       !
! then            call main_invert(qq,fhb,uu,vv,pp)  to perform inversion!
!========================================================================!

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: ld

!---------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

 !Fractional y grid values: 
fac=one/dble(ny)
do iy=0,ny
  yh1(iy)=fac*dble(iy)
  yh0(iy)=one-yh1(iy)
enddo

 !Define part of streamfunction proportional to the mean vorticity 
 !as well as beta*y:
do iy=0,ny
  yg(iy)=ymin+gly*dble(iy)
  pbar(iy)=f12*(yg(iy)-ymin)*(yg(iy)-ymax)
  bety(iy)=beta*yg(iy)
enddo

 !Define x wavenumbers:
rkx(0)=zero
do kx=1,nwxm1
  kxc=nx-kx
  rkx(kx )=hrkx(2*kx)
  rkx(kxc)=hrkx(2*kx)
enddo
rkx(nwx)=hrkx(nx)

scx=twopi/ellx
rkxmax=scx*dble(nwx)
frkx(0)=zero
do kx=1,nxm1
  wratx=rkx(kx)/rkxmax
  frkx(kx)=rkx(kx)*exp(-36.d0*wratx**36.d0)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=1,ny
  wraty=rky(ky)/rkymax
  frky(ky)=rky(ky)*exp(-36.d0*wraty**36.d0)
enddo
 
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do kx=0,nxm1
  k=nint(rkx(kx)*delki)
  kmag(kx,0)=k
  spmf(k)=spmf(k)+one
enddo
do ky=1,ny
  do kx=0,nxm1
    k=nint(sqrt(rkx(kx)**2+rky(ky)**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

 !Define inverse spectral inversion operators:
if (barot) then
   !The first mode is barotropic (kd1 = 0):
  qgop1(0,0)=zero
else
  qgop1(0,0)=-one/(kd1sq+small)
endif
qgop2(0,0)=-one/kd2sq
do kx=1,nxm1
  rksq=rkx(kx)**2
  qgop1(kx,0)=-one/(rksq+kd1sq)
  qgop2(kx,0)=-one/(rksq+kd2sq)
enddo
do ky=1,ny
  rksq=rky(ky)**2
  qgop1(0,ky)=-one/(rksq+kd1sq)
  qgop2(0,ky)=-one/(rksq+kd2sq)
enddo
do ky=1,ny
  do kx=1,nxm1
    rksq=rkx(kx)**2+rky(ky)**2
    qgop1(kx,ky)=-one/(rksq+kd1sq)
    qgop2(kx,ky)=-one/(rksq+kd2sq)
  enddo
enddo

 !Hyperbolic functions used to correct boundary conditions in inversion:
if (barot) then
   !The first mode is barotropic (kd1 = 0):
  do iy=1,nym1
    decy1(iy,0)=yh1(iy)
  enddo
  do kx=1,nxm1
    fac=rkx(kx)*elly
    div=one/(one-exp(-two*fac))
    do iy=1,nym1
      argm=fac*(one-yh1(iy))
      argp=fac*(one+yh1(iy))
      decy1(iy,kx)=(exp(-argm)-exp(-argp))*div
    enddo
  enddo
else
  do kx=0,nxm1
    fac=sqrt(rkx(kx)**2+kd1sq)*elly
    div=one/(one-exp(-two*fac))
    do iy=1,nym1
      argm=fac*(one-yh1(iy))
      argp=fac*(one+yh1(iy))
      decy1(iy,kx)=(exp(-argm)-exp(-argp))*div
    enddo
  enddo
endif

do kx=0,nxm1
  fac=sqrt(rkx(kx)**2+kd2sq)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy2(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

if (.not. barot) then
   !The first mode is *not* barotropic (kd1 .ne. 0):
  div=one/sinh(kd1*elly+small)
  ld=one/(kd1+small)
  do iy=0,ny
    pftop1(iy)= ld*cosh(kd1*(yg(iy)-ymax))*div
    pfbot1(iy)=-ld*cosh(kd1*(yg(iy)-ymin))*div
    uftop1(iy)=   -sinh(kd1*(yg(iy)-ymax))*div
    ufbot1(iy)=    sinh(kd1*(yg(iy)-ymin))*div
  enddo
  facp1=kd1sq*elly
  facu1=-kd1*elly*sinh(kd1*elly)/(cosh(kd1*elly)+small-one)
endif

div=one/sinh(kd2*elly)
ld=one/kd2
do iy=0,ny
  pftop2(iy)= ld*cosh(kd2*(yg(iy)-ymax))*div
  pfbot2(iy)=-ld*cosh(kd2*(yg(iy)-ymin))*div
  uftop2(iy)=   -sinh(kd2*(yg(iy)-ymax))*div
  ufbot2(iy)=    sinh(kd2*(yg(iy)-ymin))*div
enddo

return
end subroutine

!======================================
subroutine main_invert(qq,fhb,uu,vv,pp)

!f2py intent(in)  :: qq, fhb
!f2py intent(out) :: uu, vv, pp

 !Input:  PV field (qq) and scaled bottom topography f_0*H_b/H_1 (fhb)
 !Output: Velocity field (uu) and streamfunction field (pp)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nxm1,2),pp(0:ny,0:nxm1,2)
double precision:: uu(0:ny,0:nxm1,2),vv(0:ny,0:nxm1,2)
double precision:: fhb(0:ny,0:nxm1)
 !Local arrays:
double precision:: p1(0:ny,0:nxm1),p2(0:ny,0:nxm1),ss(0:nxm1,0:ny)
double precision:: pbot(0:nxm1),ptop(0:nxm1),cppy(nym1,0:nxm1)
double precision:: q1lin(0:ny),p1add(0:ny),u1add(0:ny)
double precision:: q2lin(0:ny),p2add(0:ny),u2add(0:ny)
double precision:: ldsq

!---------------------------------------------------------------------
! f2py hack
call init_spectral

!---------------------------------------------------------------------
 !Divide PV anomaly (q - beta*y) into modes:
do ix=0,nxm1
  do iy=0,ny
     !Add topography (f_0*H_b/H_1) to PV in lower layer:
    qq1=qq(iy,ix,1)-bety(iy)+fhb(iy,ix)
    qq2=qq(iy,ix,2)-bety(iy)
    p1(iy,ix)=vec11*qq1+vec12*qq2
    p2(iy,ix)=vec21*qq1+vec22*qq2
  enddo
enddo

 !Note: vec11 + vec12 = 1 while vec2sum = vec21 + vec22
 !Here, p1 and p2 are used temporarily

!-----------------------------------------
 !Solve for the mode 1 streamfunction, p1:

if (barot) then
   !Here, the first mode is truly barotropic (kd1 = 0); invert Laplace's 
   !operator:

   !(1) compute mean vorticity (zbar) by averaging barotropic PV anomaly:
  zbar=zero
  do ix=0,nxm1
    zbar=zbar+p1(0,ix)+p1(ny,ix)
  enddo
  zbar=f12*zbar
  do ix=0,nxm1
    do iy=1,nym1
      zbar=zbar+p1(iy,ix)
    enddo
  enddo
  zbar=zbar*dsumi

   !(2) Remove mean vorticity from p1:
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1(iy,ix)-zbar
    enddo
  enddo
 
   !(3) FFT p1 and invert to get uncorrected streamfunction p1:
  call ptospc_fc(nx,ny,p1,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=0,ny
    do kx=0,nxm1
      ss(kx,ky)=qgop1(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop_fc(nx,ny,ss,p1,xfactors,yfactors,xtrig,ytrig)

   !(4) Add part of p1 due to mean vorticity found above:
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1(iy,ix)+zbar*pbar(iy)
    enddo
  enddo

   !(5) Do a sine transform of p1 at y = ymin and ymax and obtain the
   !    interior field (cppy) that must be subtracted to give p1 = 0
   !    at y = ymin and ymax:
  do ix=0,nxm1
    pbot(ix)=p1(0,ix)
    ptop(ix)=p1(ny,ix)
  enddo
  call forfft(1,nx,pbot,xtrig,xfactors)
  call forfft(1,nx,ptop,xtrig,xfactors)

   !Define the interior semi-spectral field:
  do kx=0,nxm1
    do iy=1,nym1
      cppy(iy,kx)=pbot(kx)*decy1(ny-iy,kx)+ptop(kx)*decy1(iy,kx)
    enddo
  enddo
   !Invert using a full transform in x:
  call revfft(nym1,nx,cppy,xtrig,xfactors)

   !(6) Remove cppy to obtain the final barotropic streamfunction p1:
  do ix=0,nxm1
    p1(0, ix)=zero
    p1(ny,ix)=zero
  enddo

  do ix=0,nxm1
    do iy=1,nym1
      p1(iy,ix)=p1(iy,ix)-cppy(iy,ix)
    enddo
  enddo

   !Ensure u1bar = 0 to enforce zero global momentum:
  u1bar=zero

else
   !Here, kd1 > 0; invert the mode 1 Helmholtz operator:

   !(1) compute average mode 1 PV at each edge of the channel:
  q1bot=zero
  q1top=zero
  do ix=0,nxm1
    q1bot=q1bot+p1(0 ,ix)
    q1top=q1top+p1(ny,ix)
  enddo
  q1bot=q1bot/dble(nx)
  q1top=q1top/dble(nx)

   !(2) Subtract a linear profile from p1 before inversion:
  dq1dy=(q1top-q1bot)/elly
  fac=dq1dy*gly
  do iy=0,ny
    q1lin(iy)=q1bot+fac*dble(iy)
  enddo
  u1bar=dq1dy/(kd1sq+small)
   !u1bar: mean zonal velocity associated with this linear profile
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1(iy,ix)-q1lin(iy)
    enddo
  enddo

   !(3) FFT p1 and invert to get uncorrected streamfunction p1:
  call ptospc_fc(nx,ny,p1,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=0,ny
    do kx=0,nxm1
      ss(kx,ky)=qgop1(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop_fc(nx,ny,ss,p1,xfactors,yfactors,xtrig,ytrig)

   !(4) Do a sine transform of p1 at y = ymin and ymax and obtain the
   !    interior field (cppy) that must be subtracted to give p1 = 0
   !    at y = ymin and ymax:
  do ix=0,nxm1
    pbot(ix)=p1(0,ix)
    ptop(ix)=p1(ny,ix)
  enddo
  call forfft(1,nx,pbot,xtrig,xfactors)
  call forfft(1,nx,ptop,xtrig,xfactors)

   !Define the interior semi-spectral field:
  do kx=0,nxm1
    do iy=1,nym1
      cppy(iy,kx)=pbot(kx)*decy1(ny-iy,kx)+ptop(kx)*decy1(iy,kx)
    enddo
  enddo
   !Invert using a full transform in x:
  call revfft(nym1,nx,cppy,xtrig,xfactors)

   !(5) Remove cppy to obtain the final mode 1 streamfunction p1:
  do ix=0,nxm1
    p1(0, ix)=zero
    p1(ny,ix)=zero
  enddo

  do ix=0,nxm1
    do iy=1,nym1
      p1(iy,ix)=p1(iy,ix)-cppy(iy,ix)
    enddo
  enddo
endif

!-----------------------------------------
 !Solve for the mode 2 streamfunction, p2:

 !(1) compute average mode 2 PV at each edge of the channel:
q2bot=zero
q2top=zero
do ix=0,nxm1
  q2bot=q2bot+p2(0 ,ix)
  q2top=q2top+p2(ny,ix)
enddo
q2bot=q2bot/dble(nx)
q2top=q2top/dble(nx)

 !(2) Subtract a linear profile from p2 before inversion:
dq2dy=(q2top-q2bot)/elly
fac=dq2dy*gly
do iy=0,ny
  q2lin(iy)=q2bot+fac*dble(iy)
enddo
u2bar=dq2dy/kd2sq
 !u2bar: mean zonal velocity associated with this linear profile
do ix=0,nxm1
  do iy=0,ny
    p2(iy,ix)=p2(iy,ix)-q2lin(iy)
  enddo
enddo

 !(3) FFT p2 and invert to get uncorrected streamfunction p2:
call ptospc_fc(nx,ny,p2,ss,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    ss(kx,ky)=qgop2(kx,ky)*ss(kx,ky)
  enddo
enddo
call spctop_fc(nx,ny,ss,p2,xfactors,yfactors,xtrig,ytrig)

 !(4) Do a sine transform of p2 at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give p2 = 0
 !    at y = ymin and ymax:
do ix=0,nxm1
  pbot(ix)=p2(0,ix)
  ptop(ix)=p2(ny,ix)
enddo
call forfft(1,nx,pbot,xtrig,xfactors)
call forfft(1,nx,ptop,xtrig,xfactors)

 !Define the interior semi-spectral field:
do kx=0,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy2(ny-iy,kx)+ptop(kx)*decy2(iy,kx)
  enddo
enddo
 !Invert using a full transform in x:
call revfft(nym1,nx,cppy,xtrig,xfactors)

 !(5) Remove cppy to obtain the final mode 2 streamfunction p2:
do ix=0,nxm1
  p2(0, ix)=zero
  p2(ny,ix)=zero
enddo

do ix=0,nxm1
  do iy=1,nym1
    p2(iy,ix)=p2(iy,ix)-cppy(iy,ix)
  enddo
enddo

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
 !Re-construct layer streamfunctions by projecting modes:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix,1)=vect11*p1(iy,ix)+vect12*p2(iy,ix)
    pp(iy,ix,2)=vect21*p1(iy,ix)+vect22*p2(iy,ix)
  enddo
enddo

 !Note: these do not include the parts associated with the 
 !      linear PV profiles above; they are added back below and
 !      a further correction is made to satisfy imposed boundary
 !      conditions.

!!--------------------------------------------------------------
 !Compute velocity in each layer:
call getvel(pp,uu,vv,u1bar,u2bar)
 !Note: uu & vv are not yet corrected, though they include the
 !      contributions from the linear PV profiles above.

!!---------------------------------------------------------------
! !Correct streamfunction and velocity fields

if (.not. barot) then
   !Mode 1 is not barotropic (kd1 > 0); correct streamfunction and velocity
   !of this mode depending on choice made for mod1bc in parameters.f90:
  if (mod1bc .eq. 1) then
     !The mode 1 zonal velocity takes prescribed values of 
     !u1bot at y = ymin, and u1top at y = ymax:

     !Project zonally-averaged edge values of u onto mode 1:
    cbot=zero
    ctop=zero
    do ix=0,nxm1
      cbot=cbot+vec11*uu(0 ,ix,1)+vec12*uu(0 ,ix,2)
      ctop=ctop+vec11*uu(ny,ix,1)+vec12*uu(ny,ix,2)
    enddo
    cbot=u1bot-cbot/dble(nx)
    ctop=u1top-ctop/dble(nx)

    ldsq=one/(kd1sq+small)
    do iy=0,ny
      p1add(iy)=ctop*pfbot1(iy)+cbot*pftop1(iy)-ldsq*q1lin(iy)
      u1add(iy)=ctop*ufbot1(iy)+cbot*uftop1(iy)
    enddo
     !Above we also add back the linear part of the mode 1 streamfunction

     !Correct fields:
    do ix=0,nxm1
      do iy=0,ny
        pp(iy,ix,1)=pp(iy,ix,1)+vect11*p1add(iy)
        pp(iy,ix,2)=pp(iy,ix,2)+vect21*p1add(iy)
        uu(iy,ix,1)=uu(iy,ix,1)+vect11*u1add(iy)
        uu(iy,ix,2)=uu(iy,ix,2)+vect21*u1add(iy)
      enddo
    enddo

  else
     !Set the mean displacement of the upper free surface to zero and
     !the associated mean velocity (of mode 1) to zero:

     !Get average uncorrected streamfunction for this mode:
    call average(p1,p1avg)
     !Include contribution from linear part taken away above:
    ldsq=one/(kd1sq+small)
    p1avg=p1avg-f12*ldsq*(q1lin(0)+q1lin(ny))
     !Get average uncorrected zonal velocity for this mode (overwrite p1):
    do ix=0,nxm1
      do iy=0,ny
        p1(iy,ix)=vec11*uu(iy,ix,1)+vec12*uu(iy,ix,2)
      enddo
    enddo
    call average(p1,u1avg)

     !Enforce zero global average (mode 1) streamfunction and zonal velocity:
    aaa=facp1*p1avg
    bbb=facu1*u1avg

    cbot=f12*(bbb-aaa)
    ctop=f12*(aaa+bbb)

    do iy=0,ny
      p1add(iy)=ctop*pfbot1(iy)+cbot*pftop1(iy)-ldsq*q1lin(iy)
      u1add(iy)=ctop*ufbot1(iy)+cbot*uftop1(iy)
    enddo
     !Above we also add back the linear part of the mode 1 streamfunction

     !Correct fields:
    do ix=0,nxm1
      do iy=0,ny
        pp(iy,ix,1)=pp(iy,ix,1)+vect11*p1add(iy)
        pp(iy,ix,2)=pp(iy,ix,2)+vect21*p1add(iy)
        uu(iy,ix,1)=uu(iy,ix,1)+vect11*u1add(iy)
        uu(iy,ix,2)=uu(iy,ix,2)+vect21*u1add(iy)
      enddo
    enddo

  endif

endif

 !---------------------------------------------------
 !Next deal with mode 2 (we assume kd2 > kd1 always):

 !Project zonally-averaged edge values of u onto mode 2:
cbot=zero
ctop=zero
do ix=0,nxm1
  cbot=cbot+vec21*uu(0 ,ix,1)+vec22*uu(0 ,ix,2)
  ctop=ctop+vec21*uu(ny,ix,1)+vec22*uu(ny,ix,2)
enddo
cbot=u2bot-cbot/dble(nx)
ctop=u2top-ctop/dble(nx)

ldsq=one/kd2sq
do iy=0,ny
  p2add(iy)=ctop*pfbot2(iy)+cbot*pftop2(iy)-ldsq*q2lin(iy)
  u2add(iy)=ctop*ufbot2(iy)+cbot*uftop2(iy)
enddo
 !Above we also add back the linear part of the mode 2 streamfunction

 !Correct fields:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix,1)=pp(iy,ix,1)+vect12*p2add(iy)
    pp(iy,ix,2)=pp(iy,ix,2)+vect22*p2add(iy)
    uu(iy,ix,1)=uu(iy,ix,1)+vect12*u2add(iy)
    uu(iy,ix,2)=uu(iy,ix,2)+vect22*u2add(iy)
  enddo
enddo

 !--------------------------------------------------------------------
 !Correct zonal velocity so that x momentum h1*<u1>+h2*alpha*<u2> = 0:
call average(uu(0,0,1),u1avg)
call average(uu(0,0,2),u2avg)
umean=(h1*u1avg+alpha*h2*u2avg)*massfac
 !massfac = 1/(h1+h2*alpha) is inversely proportional to total mass

 !Subtract umean from zonal velocity component in each layer:
do lay=1,2
  do ix=0,nxm1
    do iy=0,ny
      uu(iy,ix,lay)=uu(iy,ix,lay)-umean
    enddo
  enddo
enddo

return
end subroutine

!===================================================================

subroutine getvel(pp,uu,vv,u1bar,u2bar)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -d(pp)/dy and vv = d(pp)/dx.  The mean modal velocities
! u1bar & u2bar contribute a mean flow in each layer.
! *** pp, uu & vv are all in physical space
! *** and include the domain edges.

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: pp(0:ny,0:nxm1,2),uu(0:ny,0:nxm1,2),vv(0:ny,0:nxm1,2)

 !Local arrays:
double precision:: ppi(ny  ,0:nxm1),pps(0:nxm1  ,ny)
double precision:: ppx(0:nxm1,  ny),vvi(  ny,0:nxm1)
double precision:: ppy(0:nxm1,0:ny),uui(0:ny,0:nxm1)
double precision:: ubar(2)

 !Define mean velocity in each layer associated with linear (in y) 
 !part of the zonal-average PV in each mode:
ubar(1)=vect11*u1bar+vect12*u2bar
ubar(2)=vect21*u1bar+vect22*u2bar

 !Loop over layers and compute velocity field:
do lay=1,2

   !Copy non-zero interior values of pp to ppi:
  do ix=0,nxm1
    do iy=1,nym1
      ppi(iy,ix)=pp(iy,ix,lay)
    enddo
  enddo

   !Transform ppi to spectral space:
  call ptospc_fs(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

   !Compute d(ppi)/dx = ppx spectrally:
  call xderiv_fs(nx,ny,hrkx,pps,ppx)

   !Transform ppx back to physical space as vvi:
  call spctop_fs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

   !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
  do ix=0,nxm1
    vv(0,ix,lay)=zero
    do iy=1,nym1
      vv(iy,ix,lay)=vvi(iy,ix)
    enddo
    vv(ny,ix,lay)=zero
  enddo

   !Compute d(ppi)/dy = ppy spectrally:
  call yderiv_fs(nx,ny,rky,pps,ppy)

   !Transform ppy back to physical space as uui:
  call spctop_fc(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

   !Copy -uui into uu and add mean layer velocity:
  do ix=0,nxm1
    do iy=0,ny
      uu(iy,ix,lay)=ubar(lay)-uui(iy,ix)
    enddo
  enddo

enddo
 !ends loop over layers

return
end subroutine

!===================================================================

subroutine spec1d_fc(rvar,spec)
! Computes the 1d spectrum of a field rvar which is
! periodic in x and represented by a cosine series in y
! and returns the result in spec.
!   *** Warning: rvar is modified by this routine.***

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,0:nxm1),spec(0:max(nx,ny))
 !Local array:
double precision:: wks(0:nxm1,0:ny)

 !Transform rvar to spectral space:
call ptospc_fc(nx,ny,rvar,wks,xfactors,yfactors,xtrig,ytrig)

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+f14*wks(0,0)**2

 !y-independent mode:
do kx=1,nxm1
  k=kmag(kx,0)
  spec(k)=spec(k)+f12*wks(kx,0)**2
enddo

 !x-independent mode:
do ky=1,ny
  k=kmag(0,ky)
  spec(k)=spec(k)+f12*wks(0,ky)**2
enddo

 !All other modes:
do ky=1,ny
  do kx=1,nxm1
    k=kmag(kx,ky)
    spec(k)=spec(k)+wks(kx,ky)**2
  enddo
enddo

return
end subroutine

!==================================================

end module     
