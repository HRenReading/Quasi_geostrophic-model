Module Sizes

    integer, parameter  ::  nxx =  257                ! zonal (x) dimension of the system
    integer, parameter  ::  nyy =  129                ! medridional (y) dimension of the system
    integer, parameter  ::  nl  = 2                  ! number of vertical layers (has to be 2 for this code)
    integer, parameter  ::  ntotal  = nxx*nyy*nl     ! total number of grid points
    integer, parameter  ::  nx =  288                 ! zonal (x) dimension of random field (has to be even) 
    integer, parameter  ::  ny =  144                 ! medridional (y) dimension of random field (has to be even)
    real, parameter     ::  dx = 5.e3                ! grid distance in m
    integer, parameter  ::  ntimesDA = 10             ! number of DA steps
    integer, parameter  ::  deltaDA = 24             ! time between DA steps (in model time steps)
    real, parameter     ::  dt =  3600 ! 900              ! time step in seconds
    real, parameter     ::  diff = 1.e2              ! diffusion coefficient in m^2/s in laplacian for vorticity. 
    real, parameter     ::  lengthscale = 50.         ! decorrelation lengthscale for random fields in grid points!
    real, parameter     ::  dist = 100               ! localization radius squared
    real, parameter     ::  step = 2.e-17                  ! increment step size
    integer, parameter  ::  nobs = 50*24                 ! numer of observations at each observation time
    integer, parameter  ::  nens = 500                    ! numer of ensemble members
    real, parameter     ::  sigy = 5.e3                  ! observation error standard deviation
    real, parameter     ::  Binit =  15.
    logical, parameter  ::  nrestart = 0                ! 0 is no restart 1 is restart from file 'stream000.'//resaaa  (see below)
    real, parameter     ::  w_r = 1.0
    character(5), parameter::  resaaa = '02000' 

!   arrays used by many subroutines
    real :: f(nyy),H1,H2,F1,F2,lambda
    real :: sh(0:16)
    real :: c1
    real :: cbt,sumpsibt,sumpsiboundarybt
    real :: cbc1,cbc2,n1old,n2old,n3old,a1old,a2old
    real :: xr(nx,ny),qo(nl),alpha

End Module Sizes



