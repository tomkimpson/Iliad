module derivatives

use parameters
use constants
use metric
implicit none

public geodesic, derivs

private

contains


subroutine geodesic()

 
print *, 'geodesic'
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
