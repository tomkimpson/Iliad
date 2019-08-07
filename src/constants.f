module constants

use parameters

implicit none


!Universal constants
real(kind=dp), parameter :: Newton_g = 6.67408d-11 
real(kind=dp), parameter :: Msolar = 1.989d30 
real(kind=dp), parameter :: mu = Newton_g*MBH*Msolar
real(kind=dp), parameter :: light_c = 3.0d8
real(kind=dp), parameter :: convert_m = light_c**2/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_s = light_c**3/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_spin= light_c/(Newton_g*(MBH*Msolar)**2.0_dp) !Multiply by this to go TO Natural units


!MPD
real(kind=dp), parameter :: KeplerianPeriodSeconds = KeplerianPeriod*365.0_dp*24.0_dp*3600.0_dp
real(kind=dp), parameter :: SM3 =(mu*KeplerianPeriodSeconds**2.0_dp)/(4.0_dp * PI**2.0_dp)
real(kind=dp), parameter :: semi_major = convert_m*SM3**(1.0_dp/3.0_dp)
real(kind=dp), parameter :: semi_latus = semi_major * (1.0_dp - eccentricity**2.0_dp)
real(kind=dp), parameter :: ra = semi_latus/(1.0_dp - eccentricity)
real(kind=dp), parameter :: rp = semi_latus/(1.0_dp + eccentricity)
real(kind=dp), parameter :: theta_min = (90.0_dp - iota) * PI/180.0_dp !Minimum latitude reached in radians
real(kind=dp), parameter :: zMinus = cos(theta_min)**2.0_dp
real(kind=dp), parameter :: m0 = MPSR/MBH
real(kind=dp), parameter :: inertia = 0.40_dp*(MPSR*Msolar)*(RPSR*1.0d3)**2.0_dp !SI units
real(kind=dp), parameter :: s0 = convert_spin*2.0_dp*PI*inertia/p0 !magnitude of spin spatial vector in natural units
integer(kind=dp), parameter :: entries = 12 !Number of differetnai eqns 4x(position,spin,momentum)


!Some globally defined MPD variables
real(kind=dp) :: m_sq, s_sq ! mass sqaures and s squared from initial condiitons module







end module constants
