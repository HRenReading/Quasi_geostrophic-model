Program Smooth

use Sizes

IMPLICIT NONE
 
integer time_to_start,time_to_stop,n,ntime
integer time, status,member
real  :: psi(nxx,nyy,nl,nens),psi0(nxx,nyy,nl),psit(nxx,nyy,nl),y(nobs),psiRandom(nxx,nyy,nl)
real  :: phi_g(nens),phi_r
logical new
character*9 aaa



open(77,file='output',form='formatted')

psi = 0.
call InitiateModel(psi0)
psit = psi0
 


open(10,file='stream000.02000',form='unformatted')
read(10)psit
close(10)

open(10,file='stream999.00000',form = 'unformatted')
write(10)psit
close(10)


time_to_start = 1
time_to_stop  = 10

! Generate true run and generate observations
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


call InitiateEnsemble(psi)
! Generate ensemble runs
open(10,file='phi_guess',form = 'unformatted')
read(10)phi_g
close(10)
do n = 1,nens
  !write(aaa,'(i3.3, A, i5.5)')n,'.00000'
  !open(10,file='stream'//aaa,form = 'unformatted')
  !read(10)psit
  !close(10)
  write(*,*) ' ensemble number = ', n
  do time = time_to_start, time_to_stop
    !write(*,*) ' time = ', time
    if(mod(time-1,deltaDA) .eq. 0) then
      new =.True.
    else
      new = .False.
    endif
    call IntegrateModel (psi(:,:,:,n),n,time,new,phi_g(n))
  enddo
enddo


print *, " Program is finished."

End Program Smooth

!!!!!!!!!!!!!!!!!!!!!!
subroutine InitiateEnsemble(psi)

use Sizes

IMPLICIT NONE

integer n,num,time
real psit(nxx,nyy,nl),psi(nxx,nyy,nl,nens),psiRandom(nxx,nyy,nl),phi_g(nens)
character*9 aaa
logical new

open(10,file='phi_guess',form='unformatted')
read(10)phi_g
close(10)
iloop: do n=1,nens
  open(10,file='stream000.01980',form='unformatted')
  read(10)psit
  close(10)
  !call GenerateRandomField(psiRandom)
  !psit = psit + Binit * psiRandom
  jloop: do time=1,20
    call IntegrateModel (psit(:,:,:),909,time,new,phi_g(n))
  end do jloop
  psi(:,:,:,n) = psit(:,:,:)
  write(aaa,'(i3.3, A, i5.5)')n,'.00000'
  open(10,file='stream'//aaa,form = 'unformatted')
  write(10)psi(:,:,:,n)
  close(10)
  write(*,*) ' ensemble initial = ', n
end do iloop

end subroutine InitiateEnsemble

