Program Smooth

use Sizes

IMPLICIT NONE
 
integer time_to_start,time_to_stop,n,ntime
integer time, status,member
real  :: psi(nxx,nyy,nl),psi0(nxx,nyy,nl),psit(nxx,nyy,nl),y(nobs),psiRandom(nxx,nyy,nl)
real  :: phi_g(nens),phi_r
logical new
character*3 aaa


open(77,file='output',form='formatted')

psi = 0.
call InitiateModel(psi0)
psit = psi0

time_to_start = 1
time_to_stop  = 10

! Generate true run
phi_r = exp(-1/w_r)
open(10,file='ini_r',form='unformatted')
read(10)psit
close(10)

phi_r = exp(-1/w_r)
do time = time_to_start, time_to_stop
  !write(*,*) ' time = ', time
  if(mod(time-1,deltaDA) .eq. 0) then
    new =.True.
  else
    new = .False.
  endif
  call IntegrateModel (psit(:,:,:),999,time,new,phi_r)
enddo

! Generate ensemble runs
open(10,file='phi_guess',form = 'unformatted')
read(10)phi_g
close(10)
do n = 1,nens
  write(aaa,'(i3.3, A, i5.5)')n
  open(10,file='ini'//aaa,form = 'unformatted')
  read(10)psi
  close(10)
  write(*,*) ' ensemble number = ', n
  do time = time_to_start, time_to_stop
    !write(*,*) ' time = ', time
    if(mod(time-1,deltaDA) .eq. 0) then
      new =.True.
    else
      new = .False.
    endif
    call IntegrateModel (psi(:,:,:),n,time,new,phi_g(n))
  enddo
enddo


print *, " Program is finished."

End Program Smooth

!!!!!!!!!!!!!!!!!!!!!!
subroutine measure(psi,y)

use Sizes

IMPLICIT NONE

real  :: psi(nxx,nyy,nl)
real  :: y(nobs)
integer :: i,j,k, ii

ii=0
k = 1
do i =1,50
  do j=1,24
    ii = ii+1
    y(ii) = psi(5+(i-1)*5,5+(j-1)*5,k)*psi(5+(i-1)*5,5+(j-1)*5,k)
  enddo
enddo

end subroutine measure

!!!!!!!!!!!!!!!!!!!!!!
subroutine HTRmin1innovation(psi,y,ynn,xynn)

use Sizes

IMPLICIT NONE

real psi(nxx,nyy,nl)
real  :: y(nobs),ynn(nobs),xynn(ntotal)
integer :: i,j,ix,iy,ik,n,itest,ii

xynn = 0.
itest = 0
ik = 1
do i =1,50
  ix = 5+(i-1)*5
  do j=1,24
    iy = 5+(j-1)*5
    itest = itest+1
    ii = ix + (iy-1)*nxx 
    xynn(ii) = 2*psi(ix,iy,ik)*(y(itest)-ynn(itest))/(sigy*sigy)
  enddo
enddo

end subroutine HTRmin1innovation
