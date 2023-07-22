!Subroutine GenerateRandomField 
Subroutine GenerateRandomField (psiRandom)

! Generate a correlated random field specified by the filter {xr,xi}, 
! the vertical correlation coefficient "alpha" and the amplitudes qo.

!use NagInterface        ! for UniformRandomNumbers and ComplexFFT
use Sizes               ! nxx, nyy, nl,nx,ny(somewhat larger than nxx and nyy), array xr 

IMPLICIT NONE

real psiRandom(nxx,nyy,nl)
!real, dimension(:,:,:), intent(OUT)  :: psiRandom
real, dimension(:,:), Allocatable   :: re, im, phi
real, dimension(:), allocatable     :: trigm, trign, work
!real trigm(2*nx),trign(2*ny),work(2*nx*ny)
integer                             :: ifail

integer      :: status

integer      :: i,j,k
real         :: twopi

allocate( re(nx,ny), Stat=status)      
if (status /= 0) stop "Allocation of re in randomm failed"
allocate( im(nx,ny), Stat=status)
if (status /= 0) stop "Allocation of im in randomm failed"
allocate( phi(nx,ny), Stat=status)
if (status /= 0) stop "Allocation of phi in randomm failed"

twopi  = 8. * atan(1.)
 
do k = 1, nl

  ! Fill the array phi with ion, niform randoms drawn from [0, 2 Pi)
  call random_number(phi)
  phi = phi * twopi
  !  call UniformRandomNumbers(0., twopi, phi)

! Construct phi (filter) from random numbers with the correct symmetry
  phi(1,1) = 0.

  do i = 2, nx
    do j = 2, ny/2+1
      phi(i,j)=-phi(nx-i+2,ny-j+2)
    enddo
  enddo

  do j = 2, ny/2+1
    phi(1,j)=-phi(1,ny-j+2)
  enddo

  do i = 2, nx/2+1
    phi(i,1)=-phi(nx-i+2,1)
  enddo
 
! Construct re and im
  re =  xr * cos (phi)
  im =  xr * sin (phi)
  allocate (work(2 * nx * ny), Stat = status)
  if (status /= 0) stop "Allocation of work in ComplexFFT2D failed"
  allocate (trigm(2 * nx),     Stat = status)
  if (status /= 0) stop "Allocation of trigm in ComplexFFT2D failed"
  allocate (trign(2 * ny),     Stat = status)
  if (status /= 0) stop "Allocation of trign in ComplexFFT2D failed"

  ifail = 0
  call c06fuf (nx, ny, re, im, 'i', trigm, trign, work, ifail)

  deallocate (trigm, trign, work)

! Drop additional rows/columns to reduce cyclic effects in array re.
  psiRandom(:, :, k) = re(1:nxx, 1:nyy)

enddo

deallocate (re, im, phi)

! Impose boundary condition and (re)normalize
call RestrictRandomField(psiRandom)
 
! Impose the vertical correlation
do k = 2, nl
  psiRandom(:,:,k) = sqrt(1. - alpha*alpha) * psiRandom(:, :, k) +  &
                     alpha                  * psiRandom(:, :, k-1)
enddo

! Adjust amplitudes per layer as prescribed by the input parameters qo 
do k = 1, nl
  psiRandom(:,:,k) = psiRandom(:, :, k) * qo(k)
enddo

End Subroutine GenerateRandomField


Subroutine MakeRandomFilter(r)
! Prepare the "random field generator" such that it easily creates correlated 
! field with correlation approximately equal to exp(-r * x * x), where x is 
! the distance measured in gridpoints

!use NagInterface           ! For ComplexFFT.
use Sizes                  ! sizes  nx and ny.
 
IMPLICIT NONE

real                              :: r

integer                           :: i, j, status
real                              :: x1, x2, y1, y2
real, dimension(:,:), Allocatable :: xi
real, dimension(:), allocatable   :: trigm, trign, work
integer                           :: ifail

allocate (xi(nx,ny), Stat = status)
if (status /= 0) stop "Allocation of xi in random1 failed"

do i = 1, nx
  do j = 1, ny
    x1 = i-1
    y1 = j-1
    x2 = nx - x1
    y2 = ny - y1
    xr(i,j) = x1 * y1 * exp(-r *(x2*x2+y2*y2)) + &
              x1 * y2 * exp(-r *(x2*x2+y1*y1)) + &
              x2 * y1 * exp(-r *(x1*x1+y2*y2)) + &
              x2 * y2 * exp(-r *(x1*x1+y1*y1))
  enddo
enddo

xr = xr / (nx * ny)
xi = 0. 

allocate (work(2 * nx * ny), Stat = status)
if (status /= 0) stop "Allocation of work in ComplexFFT2D failed"
allocate (trigm(2 * nx),     Stat = status)
if (status /= 0) stop "Allocation of trigm in ComplexFFT2D failed"
allocate (trign(2 * ny),     Stat = status)
if (status /= 0) stop "Allocation of trign in ComplexFFT2D failed"

ifail = 0
call c06fuf (nx, ny, xr, xi, 'i', trigm, trign, work, ifail)

deallocate (trigm, trign, work)

xr = sqrt(abs(xr))

deallocate(xi)

End Subroutine MakeRandomFilter


!Subroutine RestrictRandomField 
Subroutine RestrictRandomField (psiRandom)
! Incorporate boundary effects into random field  and scale amplitude and mean.

use Sizes               ! nxx, nyy, nl

IMPLICIT NONE

real psiRandom(nxx,nyy,nl)

integer :: i, j, k, num
real    :: average, aver2,pi

pi  = 4. * atan(1.)

psiRandom(:,1,:)   = 0.
psiRandom(:,nyy,:) = 0.
do j = 2, nyy-1
  psiRandom(:,j,:) = sin((j-2)*pi/(nyy-3)) * psiRandom(:,j,:)
enddo

do k=1,nl
  num  = COUNT (psiRandom(:,:,k) .ne. 0.)
  average = SUM (psiRandom(:,:,k)) / num

! Force the average over each layers together to be zero 
  where (psiRandom(:,:,k) .ne. 0.) psiRandom(:,:,k) = psiRandom(:,:,k) - average

  aver2 = SUM (psiRandom(:,:,k)*psiRandom(:,:,k)) / num

! Force the variance of all non-vanishing numbers to be one
  where (psiRandom(:,:,k) .ne. 0.) &
    psiRandom(:,:,k) = psiRandom(:,:,k)/ sqrt(aver2)
enddo

End Subroutine RestrictRandomField

