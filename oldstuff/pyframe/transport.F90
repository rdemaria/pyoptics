!Copyright (c) 2008, Riccardo De Maria
!All rights reserved.

#define PI 3.1415926535897931
#define LENGTH 1
#define KN0L   2
#define KS0L   3
#define KN1L   4
#define KS1L   5
#define BETAR  6
#define GAMMAR 7
#define NINPUT 7
#define SPOS   8
#define XORB   9
#define PXORB  10
#define YORB   11
#define PYORB  12
#define TORB   13
#define PTORB  14
#define DXORB  15
#define DPXORB 16
#define DYORB  17
#define DPYORB 18
#define BETX   19
#define ALFX   20
#define MUX    21
#define BETY   22
#define ALFY   23
#define MUY    24
#define NINOUT 24

function colname()
character(len=200) :: colname
colname="l kn0l ks0l kn1l ks1l br gr &
& s x px y py t pt dx dpx dy dpy betx alfx mux bety alfy muy"
end function colname


subroutine idmat(mat)
implicit none
double precision,intent(inout) :: mat(6,6)
integer :: i,j

do i=1,6
  do j=1,6
    mat(i,j)=0D0
  end do
  mat(i,i)=1D0
enddo

end subroutine


subroutine quad(v,r)
! Compute the linear transfer matrix of a
! combined function magnet
implicit none

double precision,intent(in) :: v(NINPUT)
double precision,intent(out):: r(6,6)

double precision ::l,kn0l,ks0l,kn1l,ks1l,br,gr
double precision :: cx,sx,cy,sy,dx,dy,kx,ky,skx,sky,fx

l  =v(LENGTH)
kn0l=v(KN0L)
ks0l=v(KS0L)
kn1l=v(KN1L)
ks1l=v(KS1L)
br =v(BETAR)
gr =v(GAMMAR)

r=0

if (l.eq.0d0) then
  r(1,1)=1
  r(1,2)=0
  r(2,1)=-kn1l-kn0l**2
  r(2,2)=1
  r(3,3)=1
  r(3,4)=0
  r(4,3)=kn1l-ks0l**2
  r(4,4)=1
  r(1,6)=0 !to be checked
  r(2,6)=kn0l ! to be checked
  r(5,1)=-r(2,6)
  r(5,2)=-r(1,6)
  r(3,6)=0 !to be checked
  r(4,6)=ks0l ! to be checked
  r(5,3)=-r(4,6)
  r(5,4)=-r(3,6)
  r(5,5)=1
  r(5,6)=l/br**2/gr**2
  r(6,6)=1
else
  kx=(+kn1l/l+(kn0l/l)**2)
  if (abs(kx).lt.1E-10) then
    cx=1-l**2*kx/2
    sx=l-l**3*kx/6
    dx=l**2/2
    fx=l**3/6
  else
    if (kx.gt.0D0) then
      skx=sqrt(kx)
      cx=cos(skx*l)
      sx=sin(skx*l)/skx
      dx =(1-cx)/kx
      fx =(l-sx)/kx
    else
      skx=sqrt(-kx)
      cx=cosh(skx*l)
      sx=sinh(skx*l)/skx
      dx =(1-cx)/kx
      fx =(l-sx)/kx
    endif
  endif

  ky=(-kn1l/l+(ks0l/l)**2)
  if (abs(ky).lt.1E-10) then
    cy=1-l**2*ky/2
    sy=l-l**3*ky/6
    dy=l**2/2
  else
    if (ky.gt.0D0) then
      sky=sqrt(ky)
      cy=cos(sky*l)
      sy=sin(sky*l)/sky
      dy =(1-cy)/ky
    else
      sky=sqrt(-ky)
      cy=cosh(sky*l)
      sy=sinh(sky*l)/sky
      dy =(1-cy)/ky
    endif
  endif

  r(1,1)=cx
  r(1,2)=sx
  r(2,1)=-kx*sx
  r(2,2)=cx
  r(3,3)=cy
  r(3,4)=sy
  r(4,3)=-ky*sy
  r(4,4)=cy
  r(1,6)=dx*kn0l/l
  r(2,6)=sx*kn0l/l
  r(5,1)=-r(2,6)
  r(5,2)=-r(1,6)
  r(3,6)=dy*ks0l/l
  r(4,6)=sy*ks0l/l
  r(5,3)=-r(4,6)
  r(5,4)=-r(3,6)
  r(5,5)=1
  r(5,6)=l/br**2/gr**2 -(kn0l/l)**2 * fx / br**2
  r(6,6)=1

endif

end subroutine quad


subroutine track(v,m)
implicit none

double precision,intent(inout) :: v(NINOUT,m)
integer,intent(in) :: m

double precision:: z(6)
double precision:: r(6,6),vt(NINPUT)
integer :: i


do i=1,m-1
  vt=v(1:NINPUT,i)
!  write(*,*) vt
  call quad(vt,r)
!  write(*,'(6e12.4)') r
!  write(*,*)

  ! Track s
  v(SPOS,i+1)=v(SPOS,i)+v(LENGTH,i)

  ! Track orbit
  v(XORB:PTORB,i+1)=matmul(r,v(XORB:PTORB,i))
!  ! Track dispersion
  z(1:4)=v(DXORB:DPYORB,i)
  z(5)=0D0
  z(6)=1D0
  z=matmul(r,z)
  v(DXORB:DPYORB,i+1)=z(1:4)
!  ! Track Twiss
  call trackbeta(r(1:2,1:2),v(BETX,i  ),v(ALFX,i  ),v(MUX,i  ), &
                            v(BETX,i+1),v(ALFX,i+1),v(MUX,i+1) )
  call trackbeta(r(3:4,3:4),v(BETY,i  ),v(ALFY,i  ),v(MUY,i  ), &
                            v(BETY,i+1),v(ALFY,i+1),v(MUY,i+1) )
end do

end subroutine track


subroutine tmatrix(v,m,mat)
implicit none

double precision,intent(in) :: v(NINOUT,m)
integer,intent(in) :: m

double precision,intent(out) :: mat(6,6)
double precision:: r(6,6),vt(NINPUT)
integer :: i


call idmat(mat,r)
do i=1,m-1
  vt=v(1:NINPUT,i)
!  write(*,*) vt
  call quad(vt,r)
!  write(*,'(6e12.4)') r
!  write(*,*)
  mat=matmul(mat,r)
end do

end subroutine tmatrix




subroutine trackbeta(r,bet1,alf1,mu1,bet2,alf2,mu2)
implicit none
double precision,intent(in) :: r(2,2)
double precision,intent(in) :: bet1,alf1,mu1
double precision,intent(out) :: bet2,alf2,mu2

double precision:: tmp1,tmp2,r11,r12,r21,r22

r11=r(1,1)
r12=r(1,2)
r21=r(2,1)
r22=r(2,2)

tmp1  =( r11*bet1 - r12*alf1 )
tmp2  =( r21*bet1 - r22*alf1 )
bet2 =( tmp1**2 + r12**2 ) / bet1
alf2=-( tmp1*tmp2 + r12*r22 ) / bet1
mu2=mu1 + atan2( r12,tmp1 ) / 8/atan(1D0)

end subroutine trackbeta

subroutine ptrack(v,m)
implicit none

double precision,intent(in) :: v(NINOUT,m)
integer,intent(in) :: m

double precision :: mat(6,6)

call tmatrix(v,m,mat)
call pbeta(mat(1:2,1:2),v(BETX,1),v(ALFX,1),v(MUX,1))
call pbeta(mat(3:4,3:4),v(BETY,1),v(ALFY,1),v(MUY,1))
call pdisp(mat(1:2,1:6),v(DXORB,1),v(DPXORB,1))
call pdisp(mat(3:4,1:6),v(DYORB,1),v(DPYORB,1))
call track(v,m)

end subroutine ptrack


subroutine trace(v,m,t)
implicit none

double precision,intent(in) :: v(NINOUT,m)
integer,intent(in) :: m
double precision,intent(out) :: t(6)

double precision :: mat(6,6)
integer :: i

call tmatrix(v,m,mat)
do i=1,6
  t(i)=mat(i,i)
end do

end subroutine trace


subroutine pbeta(r,bet,alf,mu)
implicit none
double precision,intent(in) :: r(2,2)
double precision,intent(out) :: bet,alf,mu

double precision :: cmu, smu,r11,r12,r21,r22

r11=r(1,1)
r12=r(1,2)
r21=r(2,1)
r22=r(2,2)

cmu=(r11+r22)/2
smu=sign(1D0,r12) * sqrt(-r12*r21 - (r11-r22)**2 / 4.)
!write(*,*) cmu,smu
mu=.5/PI * atan(smu/cmu)
bet=r12/smu
alf=.5*(r22-r11)/smu

end subroutine pbeta

subroutine pdisp(r,z)
implicit none
double precision,intent(in) :: r(2,6)
double precision,intent(out) :: z(2)

double precision:: r16,r26,r11,r12,r21,r22,det

r11=r(1,1)-1D0
r12=r(1,2)
r21=r(2,1)
r22=r(2,2)-1D0
r16=r(1,6)
r26=r(2,6)

det=r11*r22-r12*r21
if (abs(det).gt.1E-15) then
  z(1)=(-r16*r22+r26*r12) / det
  z(2)=(r26*r11-r16*r21) / det
endif



end subroutine pdisp
