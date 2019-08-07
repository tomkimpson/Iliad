module OrbitalDynamics

use parameters
use constants
use metric

implicit none

private calculate_EQL, function_f, function_g, function_h, function_d, &
        calculate_p, calculate_s

public MPD

contains


subroutine MPD()

real(kind=dp), dimension(4) :: xi !Initial spatial coordinates t_init, r_init, theta_init, phi_init, t_init
real(kind=dp), dimension(3) :: ci !Conserved quantities
real(kind=dp), dimension(4) :: pvector,svector !4-momentum
real(kind=dp) :: E,Q,L
real(kind=dp), dimension(entries) :: y_init !the array of initial conditions to pass to the rk solver



print *, 'MPD subroutine called'

xi(1) = 0.0_dp !t_init
xi(2) = semi_major
xi(3) = PI/2.0_dp
xi(4) = 0.0_dp



!Calculate the initial E,L,Q from the Keplerian orbital parameters
call calculate_EQL(E,Q,L)
ci(1) = E
ci(2) = Q
ci(3) = L

!Calculate initial 4-momentum
call calculate_p(xi,ci,pvector)

!Calculate the initial 4-spin
call calculate_s(xi,ci,pvector,svector)


!Calculate the two conserved quantities, m and s
call mag_4(xi(2),xi(3),pvector,m_sq)
m_sq = -m_sq

call mag_4(xi(2),xi(3), svector,s_sq)

!Pass it all to an array and pass that array to the numerical integrator
y_init(5:8) = xi
y_init(5:8) = pvector
y_init(9:12) = svector



end subroutine MPD


!!!!!!!!!-----PRIVATE---!!!!!!!!!1



subroutine calculate_s(xi,ci,pvector,svector)

!Arguments
real(kind=dp), dimension(4),intent(in) :: xi
real(kind=dp), dimension(3),intent(in) :: ci
real(kind=dp), dimension(4),intent(in) :: pvector
real(kind=dp), dimension(4),intent(out) :: svector

!Other
real(kind=dp) :: r,theta
real(kind=dp), dimension(4,4) :: metric


r = xi(2)
theta = xi(3)

svector(2) = s0 * sin(stheta) * cos(sphi) !S1
svector(3) = s0*sin(stheta)*sin(sphi)/r  !S2
svector(4) = s0 *cos(stheta)/(r*sin(theta)) !S3


call calculate_covariant_metric(r,theta,metric)




svector(1) = -( &
             (metric(2,1) * pvector(1) + metric(2,2)*pvector(2) + metric(2,3)*pvector(3) + metric(2,4) * pvector(4) ) * svector(2)+&
             (metric(3,1) * pvector(1) + metric(3,2)*pvector(2) + metric(3,3)*pvector(3) + metric(3,4) * pvector(4) ) * svector(3)+&
             (metric(4,1) * pvector(1) + metric(4,2)*pvector(2) + metric(4,3)*pvector(3) + metric(4,4) * pvector(4) ) * svector(4) &
             ) / &
             (metric(1,1)*pvector(1) + metric(2,1) *pvector(2) + metric(3,1)*pvector(3) + metric(4,1) * pvector(4))





end subroutine calculate_s




subroutine calculate_p(xi,ci,pvector)
!Arguments
real(kind=dp), dimension(4),intent(in) :: xi
real(kind=dp), dimension(3),intent(in) :: ci
real(kind=dp), dimension(4),intent(out) :: pvector
!Other
real(kind=dp) :: r,theta,phi,sigma,delta,RR,PP,TT, E,Q,L

!Read it in
r = xi(2)
theta = xi(3)
phi = xi(4)
E = ci(1)
Q = ci(2)
L = ci(3)


!Do some calculations
sigma = r**2.0_dp +a**2.0_dp * cos(theta)
delta = r**2.0_dp +a**2.0_dp - 2.0_dp*r


!Remember - these only apply for Kerr
PP = E* (r**2.0_dp + a**2.0_dp) - a*L
RR = PP**2.0_dp -delta*(r**2.0_dp + Q + (L-a*E)**2.0_dp)
TT = Q - cos(theta)**2.0_dp*(a**2.0_dp * (1.0_dp - E**2.0_dp)+L**2.0_dp/sin(theta)**2.0_dp)



pvector(1) = m0 * (a*(L-a*E*sin(theta)**2.0_dp) + (r**2.0_dp + a**2.0_dp)*PP/delta)/sigma
pvector(2) = m0*sqrt(RR) / sigma
pvector(3) = -m0*sqrt(TT) / sigma
pvector(4) = m0*((L/sin(theta)**2 - a*E) + a*PP/delta)/sigma




end subroutine calculate_p












subroutine calculate_EQL(E, Q, L)
real(kind=dp) :: f1,g1,h1,d1 !f_functions used in defining the determinants
real(kind=dp) :: f2,g2,h2,d2 !f_functions used in defining the determinants
real(kind=dp) :: kappa, epsil, rho, eta, sigma !determinants
real(kind=dp) :: DD ! labels prograge or retrograde orbits

real(kind=dp) :: E1, E2, E3, E !Different parts of the energy expression
real(kind=dp) :: L1, L2, L !Different parts of the momentum expression
real(kind=dp) :: Q !Carter constant


if (a .LT. 0.0_dp) then
  DD = -1.0_dp

else if (a .GT. 0.0_dp) then
  DD = 1.0_dp
else if (A .EQ. 0.0_dp) then
  DD = 1.0_dp
  print *, 'Spin parameter a = 0 (Schwarzchild). Setting DD = +1 in initial_EQL module'
  endif




call function_f(rp,f1)
call function_g(rp,g1)
call function_h(rp,h1)
call function_d(rp,d1)


call function_f(ra,f2)
call function_g(ra,g2)
call function_h(ra,h2)
call function_d(ra,d2)


kappa = d1*h2 - d2*h1
epsil = d1*g2 - d2*g1
rho = f1*h2 - f2*h1
eta = f1*g2 - f2*g1
sigma = g1*h2 - g2*h1


!Now calculate the energy

E1 = sigma*(sigma*epsil**2.0_dp + rho*epsil*kappa - eta*kappa**2.0_dp)
E2 = kappa*rho + 2.0_dp*epsil*sigma - 2.0_dp*DD*sqrt(E1)
E3 = rho**2.0_dp + 4.0_dp*eta*sigma

E = sqrt(E2/E3)




!And the angular momentum
L1 = -g1*E/h1
L2 = g1**2.0_dp * E**2.0_dp + (f1*E**2.0_dp - d1)*h1

L = L1 + DD*sqrt(L2)/h1



!And finally the Carter constant

Q = zMinus * (a**2.0_dp * (1.0_dp - E**2.0_dp) + L**2.0_dp/(1.0_dp-zMinus))

end subroutine calculate_EQL








subroutine function_f(r,f)
real(kind=dp) :: r,f ! in , out
real(kind=dp) :: delta
real(kind=dp) :: ans1, ans2, rr
integer(kind=dp) :: pow
delta = r**2 - 2.0_dp*r + a**2.0_dp
f = r**4 + a**2.0_dp * (r*(r+2.0_dp) + zMinus*delta)


end subroutine function_f


subroutine function_g(r,g)
real(kind=dp) :: r,g
g = 2.0_dp*a*r
end subroutine function_g


subroutine function_h(r,h)
real(kind=dp) :: r,h
real(kind=dp) :: delta

delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp
h = r*(r-2.0_dp) + zMinus*delta/(1.0_dp - zMinus)

end subroutine function_h




subroutine function_d(r,d)
real(kind=dp) :: r, d
real(kind=dp) :: delta

delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp
d = (r**2.0_dp + a**2.0_dp * zMinus)*delta

end subroutine function_d


















end module OrbitalDynamics
