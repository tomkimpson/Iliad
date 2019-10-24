module parameters
        
implicit none

!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)



!Universal constants
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 
!BH parameters
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: a = 0.9980_dp !BH spin parameter


!PSR parameters

real(kind=dp), parameter :: MPSR = 1.40_dp !PSR mass is solar masses
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in km
real(kind=dp), parameter :: stheta = 0.0_dp, sphi = 0.0_dp !Initial angle of spin axis
real(kind=dp), parameter :: p0 = 1.0d-3 !pulsar spin period in seconds
real(kind=dp), parameter :: chi = PI/2.0_dp !Polar angle between spin and radiation axis

!Orbital parameters
real(kind=dp), parameter :: KeplerianPeriod = 0.10_dp !years
real(kind=dp), parameter :: eccentricity = 0.10_dp
real(kind=dp), parameter :: iota = 0.0_dp !Inclination w.r.t equatorial plane in degrees 0.60_dp
real(kind=dp), parameter :: lambda = 0.0_dp !Turn on/off spin curvature coupling
real(kind=dp), parameter :: Norbit = 3.0_dp !Number of orbits


!Plasma paramters
real(kind=dp), parameter :: N = 0.0e7_dp !plasma density normalisation



!Integration options
integer(kind=dp)  :: adaptive = 1 !turn on/off (1/0) adaptive stepsize.
!Note if off, it is important to pay attention to stepsize. Not a paramter since changed between MPD and RT
real(kind=dp), parameter :: hs = 1.0d-5 !Fixed timing resolution in seconds. only used if adaptive =0
real(kind=dp), parameter :: RepRes = 1000.0_dp !The representative resolution. What fraction of the orbit do you want for target &
!points? 1 = all of the orbit.
integer(kind=dp), parameter :: RayTracingDirection = -1 !+1 = Forward, -1 = Backward
real(kind=dp), parameter :: RObs = 1000.0_dp,ThetaObs = PI/2.0_dp, PhiObs = 0.0_dp

!I/O options
integer(kind=dp), parameter :: plot_MPD = 1 !turn on/off (1/0) numerical accuracy evaluation
integer(kind=dp), parameter :: plot_RT = 1 !turn on/off (1/0) numerical accuracy evaluation


end module parameters
