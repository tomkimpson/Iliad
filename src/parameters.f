module parameters
        
implicit none

!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)



!Universal constants
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 
!BH parameters
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: a = 0.60_dp !BH spin parameter


!PSR parameters

real(kind=dp), parameter :: MPSR = 1.40_dp !PSR mass is solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in km
real(kind=dp), parameter :: stheta = PI/2.0_dp, sphi = 0.0_dp !Angle of spin axis
real(kind=dp), parameter :: p0 = 1.0d-3 !pulsar spin period in seconds


!Orbital parameters
real(kind=dp), parameter :: KeplerianPeriod = 0.10_dp !years
real(kind=dp), parameter :: eccentricity = 0.60_dp
real(kind=dp), parameter :: iota = 10.0_dp !Inclination w.r.t equatorial plane in degrees 0.60_dp




!Plasma paramters
real(kind=dp), parameter :: N = 3.50d7 !plasma density normalisation


!I/O options
integer(kind=dp), parameter :: plot = 1 !turn on/off (1/0) numerical accuracy evaluation



end module parameters
