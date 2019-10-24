module derivatives

use parameters
use constants
use metric
implicit none

public geodesic, derivs

private

contains


subroutine geodesic(v,constants,dv)


        !Arguments
real(kind=dp), dimension(6) :: v,dv
real(kind=dp), dimension(4),intent(in) :: constants
!Other
real(kind=dp) :: r,theta,phi,t,pr,ptheta
real(kind=dp) :: sigma,delta, SD,csc
real(kind=dp) :: Lz, kappa,B2, plasma_correction_pr, plasma_correction_pt
real(kind=dp) :: fr, fr_prime, ft, ft_prime

!Read in coordinate variables
r= v(1)
theta= v(2)
phi= v(3)
t= v(4)
pr= v(5)
ptheta= v(6)

!Get the constants
Lz = constants(1)
kappa = constants(2)
B2 = constants(3)

!Define some useful quantities


sigma = r**2 + a**2 * cos(theta)**2
delta = r**2 -2.0_dp*r +a**2
csc = 1.0_dp /sin(theta)
SD = sigma*delta

dv(1) = pr*delta/sigma
dv(2) = ptheta/sigma
dv(3) = (2.0_dp*a*r + (sigma - 2.0_dp*r)*Lz*csc**2 ) / SD
dv(4) = 1.0_dp + (2.0_dp*r*(r**2+a**2)-2.0_dp*a*r*Lz)/SD


call plasma_fr(r,fr)
call plasma_fr_deriv(r,fr_prime)

call plasma_fr(theta,ft)
call plasma_ft_deriv(theta,ft_prime)



plasma_correction_pr = -B2*(fr_prime*delta/2.0_dp + fr*(r-1.0_dp)) / SD
plasma_correction_pt = -B2*ft_prime/(2.0_dp*sigma)


dv(5) = ( &
        -kappa*(r-1.0_dp) &
        +2.0_dp*r*(r**2+a**2) &
        -2.0_dp*a*Lz &
        )/SD &
        - 2.0_dp*pr**2*(r-1.0_dp)/sigma + plasma_correction_pr





dv(6) = sin(theta)*cos(theta)*(Lz**2 * csc**4 -a**2) / sigma + plasma_correction_pt

 
end subroutine geodesic





subroutine derivs(y, dy)

!Arguments
real(kind=dp), intent(IN), dimension(entries) :: y
real(kind=dp), intent(OUT), dimension(entries) :: dy !deriatives

!Other
real(kind=dp), dimension(4,4) :: metric,metricCONTRA !covariant and contravariant metric components
real(kind=dp), dimension(4) ::  SVector, PVector
real(kind=dp), dimension(4) :: Xprime, Sprime, Pprime !derivative vectors
real(kind=dp) :: r,theta

!Read in the data
r = y(2)
theta = y(3)
PVector = y(5:8)
SVector = y(9:12)


!!Calculate the metric components
call calculate_covariant_metric(r,theta, metric)
call calculate_contravariant_metric(r,theta, metricCONTRA)


!Calculate Christoffel symbols - these are saved globally

call calculate_christoffel(r,theta)

!Calculate Riemann tensor - components are saved globally
call calculate_riemann(r, theta)

!Calculate components of antisymmetric spin tensor - cpts saved globally

call calculate_spintensor(r, theta, &
                            SVector, PVector, metricCONTRA)


!!Calculate 4-velocity
call calculate_FourVelocity(PVector, metric,Xprime)
!
!
!!Calculate 4-momentum
call calculate_FourMom(Xprime,PVector,Pprime)
!
!
!!Calculate 4-spin
call calculate_FourSpin(Xprime,PVector,SVector,Sprime)
!!


dy(1:4) = Xprime
dy(5:8) = Pprime
dy(9:12) = Sprime



end subroutine derivs




end module derivatives
