module metric


use parameters
use constants

implicit none

private

public calculate_contravariant_metric, calculate_covariant_metric, mag_4, calculate_christoffel, &
        calculate_riemann, calculate_spintensor, calculate_FourVelocity, calculate_FourMom, calculate_FourSpin, &
        plasma_fr, plasma_fr_deriv, plasma_ft, plasma_ft_deriv

contains




subroutine mag_4(r,theta,vector,mag)
!Arguments
real(kind=dp), intent(in) :: r,theta
real(kind=dp), intent(in), dimension(4) :: vector
real(kind=dp), intent(out) :: mag
!Other
real(kind=dp), dimension(4,4) :: metric

!WE ASSUME THE VECTOR IS CONTRAVARIANT


call calculate_covariant_metric(r,theta,metric)

mag = metric(1,1)*vector(1)**2.0_dp + metric(2,2)*vector(2)**2.0_dp + metric(3,3)*vector(3)**2.0_dp+metric(4,4)*vector(4)**2.0_dp&
      +2.0_dp*(metric(1,2)*vector(1)*vector(2) + & 
               metric(1,3)*vector(1)*vector(3) + &
               metric(1,4)*vector(1)*vector(4) + &
               metric(2,3)*vector(2)*vector(3) + &
               metric(2,4)*vector(2)*vector(4) + &
               metric(3,4)*vector(3)*vector(4)   &
               )





end subroutine mag_4


subroutine calculate_contravariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2.0_dp*cos(theta)**2
delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp


metric(1,1) = - ((r**2.0_dp + a**2.0_dp) + 2.0_dp*r*a**2.0_dp*sin(theta)**2.0_dp/sigma)/delta
metric(2,2) = delta/sigma
metric(3,3) = 1.0_dp/sigma
metric(4,4) = (1.0_dp-2.0_dp*r/sigma)/(delta*sin(theta)**2.0_dp)

metric(1,4) = -2.0_dp*r*a/(sigma*delta)
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_contravariant_metric


subroutine calculate_covariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2.0_dp*cos(theta)**2
delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp



metric(1,1) = -(1.0_dp-2.0_dp*r/sigma)
metric(2,2) = sigma/delta
metric(3,3) = sigma
metric(4,4) = (r**2.0_dp + a**2.0_dp + 2.0_dp*r*a**2.0_dp*sin(theta)**2.0_dp/sigma)*sin(theta)**2.0_dp


metric(1,4) = -2.0_dp*r*a*sin(theta)**2.0_dp/sigma
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_covariant_metric

SUBROUTINE calculate_christoffel(r,theta)
!Arguments
real(kind=dp), intent(in) :: r,theta
!Other
real(kind=dp) :: Sg, Dl


Sg = r**2+a**2*cos(theta)**2  
Dl = r**2 + a**2 - 2.0_dp*r

G0_00 = 0.0_dp
G0_01 = (r**2+a**2)*(r**2-a**2*cos(theta)**2)/(Sg**2*Dl)
G0_02 = -2.0_dp*r*a**2*sin(theta)*cos(theta)/Sg**2
G0_03 = 0.0_dp
G0_11 = 0.0_dp
G0_12 = 0.0_dp
G0_13 = -a*((3.0_dp*r**2-a**2)*(r**2+a**2)-a**2*(r**2-a**2)*sin(theta)**2)        &
*sin(theta)**2/(Sg**2*Dl)
G0_22 = 0.0_dp
G0_23 = 2.0_dp*r*a**3*sin(theta)**3*cos(theta)/Sg**2
G0_33 = 0.0_dp

G1_00 = Dl*(r**2-a**2*cos(theta)**2)/Sg**3
G1_01 = 0.0_dp
G1_02 = 0.0_dp
G1_03 = -a*Dl*(r**2-a**2*cos(theta)**2)*sin(theta)**2/Sg**3
G1_11 = (-(r**2-a**2)+a**2*sin(theta)**2*(r-1.0_dp))/(Sg*Dl)
G1_12 = -a**2*sin(theta)*cos(theta)/Sg
G1_13 = 0.0_dp
G1_22 = -r*Dl/Sg
G1_23 = 0.0_dp
G1_33 = -(r*(a**2+r**2)**2-a**2*sin(theta)**2*(-a**2*sin(theta)**2*(r-1.0_dp)            &
-1.0_dp*a**2+r*(2.0_dp*(r**2+a**2)+r)))*Dl*sin(theta)**2/Sg**3

G2_00 = -2.0_dp*r*a**2*sin(theta)*cos(theta)/Sg**3
G2_01 = 0.0_dp
G2_02 = 0.0_dp
G2_03 = 2.0_dp*r*a*(r**2+a**2)*sin(theta)*cos(theta)/Sg**3
G2_11 = a**2*sin(theta)*cos(theta)/(Sg*Dl)
G2_12 = r/Sg
G2_13 = 0.0_dp
G2_22 = -a**2*sin(theta)*cos(theta)/Sg
G2_23 = 0.0_dp
G2_33 = -((r**2+a**2)**3-a**2*Dl*(2.0_dp*(r**2+a**2)-a**2*sin(theta)**2)*sin(theta)**2) &
*sin(theta)*cos(theta)/Sg**3

G3_00 = 0.0_dp
G3_01 = a*(r**2-a**2*cos(theta)**2)/(Sg**2*Dl)
G3_02 = -2.0_dp*r*a*(cos(theta)/sin(theta))/Sg**2
G3_03 = 0.0_dp
G3_11 = 0.0_dp
G3_12 = 0.0_dp
G3_13 = (r*Dl*(r**2+a**2)+a**4*(r-1.0_dp)*sin(theta)**4-a**2*(a**2+r**2)*(2.0_dp*r-1.0_dp)*sin(theta)**2) &
/(Sg**2*Dl)
G3_22 = 0.0_dp
G3_23 = ((r**2+a**2)**2+sin(theta)**4*a**4-2.0_dp*a**2*(r**2-1.0_dp*r+a**2)*sin(theta)**2)*(cos(theta)/sin(theta))/Sg**2 
G3_33 = 0.0_dp





END SUBROUTINE calculate_christoffel


SUBROUTINE calculate_riemann(r,theta)
!Argument
real(kind = dp) :: r,theta
!Other
real(kind=dp) :: M,Dl,Sg


M = 1.0_dp
Sg = r**2+a**2*cos(theta)**2  
Dl = r**2 + a**2 - 2.0_dp*M*r

R_0101 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*Dl+a**2*sin(theta)**2)/(Sg**3*Dl)


R_0102 = 3.0_dp*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**3


R_0103 = 0.0_dp


R_0112 = 0.0_dp


R_0113 = -M*r*a*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2  &
/(Sg**3*Dl)


R_0123 = M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)           &
*sin(theta)*cos(theta)/Sg**3


R_0202 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(Dl+2.0_dp*a**2*sin(theta)**2)/Sg**3


R_0203 = 0.0_dp


R_0212 = 0.0_dp


R_0213 = M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)             &
*sin(theta)*cos(theta)/Sg**3


R_0223 = M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**3


R_0303 = M*r*Dl*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**3


R_0312 = -M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**2


R_0313 = 0.0_dp


R_0323 = 0.0_dp


R_1212 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg*Dl)


R_1213 = 0.0_dp


R_1223 = 0.0_dp


R_1313 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*((r**2+a**2)**2+2.0_dp*a**2*Dl*sin(theta)**2)  &
/(Sg**3*Dl)


R_1323 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)**3*cos(theta)/Sg**3


R_2323 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*Dl*sin(theta)**2)  &
/Sg**3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

R_0201 = R_0102


R_0301 = R_0103


R_1201 = R_0112


R_1301 = R_0113


R_2301 = R_0123


R_0302 = R_0203


R_1202 = R_0212


R_1302 = R_0213


R_2302 = R_0223


R_1203 = R_0312


R_1303 = R_0313


R_2303 = R_0323


R_1312 = R_1213


R_2312 = R_1223


R_2313 = R_1323


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R0_001 = 0.0_dp


R0_002 = 0.0_dp


R0_003 = 2.0_dp*sin(theta)**2*M**2*a*r**2*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**4


R0_012 = -2.0_dp*(3.0_dp*r**2-a**2*cos(theta)**2)*r*M**2*cos(theta)*a**2*sin(theta)/(Sg**3*Dl)


R0_013 = 0.0_dp


R0_023 = 0.0_dp


R0_101 = M*r*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg**3*Dl)


R0_102 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2*M*r)*M*a**2*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_103 = 0.0_dp


R0_112 = 0.0_dp


R0_113 = 3.0_dp*sin(theta)**2*a*M*r*(r**2+a**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/(Sg**3*Dl)


R0_123 = -cos(theta)*sin(theta)*M*a*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*sin(theta)**2*Dl)/(Sg**3*Dl)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R0_201 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4*M*r)*a**2*sin(theta)*cos(theta)*M/(Sg**3*Dl)


R0_202 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(r**2+a**2+2*a**2*sin(theta)**2)/Sg**3


R0_203 = 0.0_dp


R0_212 = 0.0_dp


R0_213 = -(3.0_dp*r**2-a**2*cos(theta)**2)*((a**2+r**2)**2+2.0_dp*a**2*sin(theta)**2*Dl)*M*a*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_223 = -3.0_dp*sin(theta)**2*M*r*a*(a**2+r**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R0_301 = 0.0_dp


R0_302 = 0.0_dp


R0_303 = -sin(theta)**2*M*r*(r**2+3*a**2*sin(theta)**2-3*a**2)*((a**2+r**2)**2-a**2*sin(theta)**2*Dl)/Sg**4


R0_312 = (3.0_dp*r**2-a**2*cos(theta)**2)*((a**2+r**2)**2-a**2*sin(theta)**2*Dl)*M*a*sin(theta)*cos(theta)/(Sg**3*Dl)


R0_313 = 0.0_dp


R0_323 = 0.0_dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R1_001 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(a**2*sin(theta)**2+2.0_dp*Dl)/Sg**4


R1_002 = -3.0_dp*Dl*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**4


R1_003 = 0.0_dp


R1_012 = 0.0_dp


R1_013 = M*r*a*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**4


R1_023 = -a*M*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl*cos(theta)*sin(theta)/Sg**4


R1_101 = 0.0_dp


R1_102 = 0.0_dp


R1_103 = 0.0_dp


R1_112 = 0.0_dp


R1_113 = 0.0_dp


R1_123 = 0.0_dp


R1_201 = 0.0_dp


R1_202 = 0.0_dp


R1_203 = -cos(theta)*sin(theta)*a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl/Sg**3


R1_212 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**2


R1_213 = 0.0_dp


R1_223 = 0.0_dp


R1_301 = -sin(theta)**2*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)/Sg**4


R1_302 = a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)*Dl*cos(theta)*sin(theta)/Sg**4


R1_303 = 0.0_dp


R1_312 = 0.0_dp


R1_313 = -sin(theta)**2*M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*((r**2+a**2)**2+2.0_dp*a**2*Dl*sin(theta)**2)/Sg**4


R1_323 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*Dl*cos(theta)*sin(theta)**3/Sg**4 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R2_001 = -3.0_dp*M*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*sin(theta)*cos(theta)/Sg**4


R2_002 = -M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*a**2*sin(theta)**2+Dl)/Sg**4


R2_003 = 0.0_dp


R2_012 = 0.0_dp


R2_013 = -a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)*cos(theta)*sin(theta)/Sg**4


R2_023 = -M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**4


R2_101 = 0.0_dp


R2_102 = 0.0_dp


R2_103 = a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/Sg**3


R2_112 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)/(Dl*Sg**2)


R2_113 = 0.0_dp


R2_123 = 0.0_dp


R2_201 = 0.0_dp


R2_202 = 0.0_dp


R2_203 = 0.0_dp


R2_212 = 0.0_dp


R2_213 = 0.0_dp


R2_223 = 0.0_dp


R2_301 = a*M*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/Sg**4


R2_302 = M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*sin(theta)**2/Sg**4


R2_303 = 0.0_dp


R2_312 = 0.0_dp


R2_313 = 3.0_dp*M*a**2*(r**2+a**2)*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)**3/Sg**4


R2_323 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(2.0_dp*(r**2+a**2)**2+a**2*Dl*sin(theta)**2)*sin(theta)**2/Sg**4 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


R3_001 = 0.0_dp


R3_002 = 0.0_dp


R3_003 = 0.0_dp


R3_012 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(a**2*sin(theta)**2-Dl)*M*a*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_013 = 0.0_dp


R3_023 = 0.0_dp


R3_101 = 3.0_dp*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)/(Dl*Sg**3)


R3_102 = -a*M*(3.0_dp*r**2-a**2*cos(theta)**2)*(2.0_dp*a**2*sin(theta)**2+Dl)*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_103 = 0.0_dp


R3_112 = 0.0_dp


R3_113 = M*r*(r**2-3.0_dp*a**2*cos(theta)**2)*(r**2+a**2+2.0_dp*a**2*sin(theta)**2)/(Dl*Sg**3)


R3_123 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-2.0_dp*M*r)*M*a**2*sin(theta)*cos(theta)/(Dl*Sg**3)


R3_201 = -(3.0_dp*r**2-a**2*cos(theta)**2)*(a**2*sin(theta)**2+2.0_dp*Dl)*M*a*(cos(theta)/sin(theta))/(Dl*Sg**3)


R3_202 = -3.0_dp*M*r*a*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R3_203 = 0.0_dp


R3_212 = 0.0_dp


R3_213 = -M*a**2*(3*r**2-a**2*cos(theta)**2)*(3.0_dp*(r**2+a**2)-4.0_dp*M*r)*sin(theta)*cos(theta)/(Dl*Sg**3)


R3_223 = -M*r*(2.0_dp*(r**2+a**2)+a**2*sin(theta)**2)*(r**2-3.0_dp*a**2*cos(theta)**2)/Sg**3


R3_301 = 0.0_dp


R3_302 = 0.0_dp


R3_303 = -2.0_dp*M**2*a*r**2*(r**2-3.0_dp*a**2*cos(theta)**2)*sin(theta)**2/Sg**4


R3_312 = 2.0_dp*M**2*r*a**2*(3.0_dp*r**2-a**2*cos(theta)**2)*cos(theta)*sin(theta)/(Dl*Sg**3)


R3_313 = 0.0_dp


R3_323 = 0.0_dp

END SUBROUTINE calculate_riemann




subroutine calculate_spintensor(r, theta, &
                                  SVector, PVector, metric)

!Arguments
real(kind=dp) :: r,theta
!Other
real(kind=dp), dimension(4,4) :: metric !This is actually the contravariant metric
real(kind=dp), dimension(4) ::  SVector, PVector
real(kind=dp) ::  S0,S1,S2,S3,& !spin components
                  P0,P1,P2,P3,& !p components
                  H00,H01,H02,H03,& !metric components
                  H10,H11,H12,H13,&
                  H20,H21,H22,H23,&
                  H30,H31, H32,H33

real(kind=dp) :: Sg, eta0123 !sigma and eta defined



!Read in component values
S0 = SVector(1)
S1 = SVector(2)
S2 = SVector(3)
S3 = SVector(4)


P0 = PVector(1)
P1 = PVector(2)
P2 = PVector(3)
P3 = PVector(4)


H00 = metric(1,1)
H01 = metric(1,2)
H02 = metric(1,3)
H03 = metric(1,4)


H10 = metric(2,1)
H11 = metric(2,2)
H12 = metric(2,3)
H13 = metric(2,4)


H20 = metric(3,1)
H21 = metric(3,2)
H22 = metric(3,3)
H23 = metric(3,4)


H30 = metric(4,1)
H31 = metric(4,2)
H32 = metric(4,3)
H33 = metric(4,4)


!And finally calculate
Sg = r**2+a**2*cos(theta)**2  
eta0123 = Sg*sin(theta)


S_01 = eta0123/m0*((H00*H11-H01*H01)*(P2*S3-P3*S2)                    &
+ (H00*H12-H01*H02)*(P3*S1-P1*S3) + (H00*H13-H01*H03)*(P1*S2-P2*S1)   & 
+ (H10*H12-H11*H02)*(P0*S3-P3*S0) + (H10*H13-H11*H03)*(P2*S0-P0*S2)   &
+ (H20*H13-H21*H03)*(P0*S1-P1*S0))

S_02 = eta0123/m0*((H00*H21-H02*H01)*(P2*S3-P3*S2)                    &
+ (H00*H22-H02*H02)*(P3*S1-P1*S3) + (H00*H23-H02*H03)*(P1*S2-P2*S1)   & 
+ (H10*H22-H12*H02)*(P0*S3-P3*S0) + (H10*H23-H12*H03)*(P2*S0-P0*S2)   &
+ (H20*H23-H22*H03)*(P0*S1-P1*S0))

S_03 = eta0123/m0*((H00*H31-H03*H01)*(P2*S3-P3*S2)                    &
+ (H00*H32-H03*H02)*(P3*S1-P1*S3) + (H00*H33-H03*H03)*(P1*S2-P2*S1)   & 
+ (H10*H32-H13*H02)*(P0*S3-P3*S0) + (H10*H33-H13*H03)*(P2*S0-P0*S2)   &
+ (H20*H33-H23*H03)*(P0*S1-P1*S0))

S_12 = eta0123/m0*((H01*H21-H02*H11)*(P2*S3-P3*S2)                    &
+ (H01*H22-H02*H12)*(P3*S1-P1*S3) + (H01*H23-H02*H13)*(P1*S2-P2*S1)   & 
+ (H11*H22-H12*H12)*(P0*S3-P3*S0) + (H11*H23-H12*H13)*(P2*S0-P0*S2)   & 
+ (H21*H23-H22*H13)*(P0*S1-P1*S0))

S_13 = eta0123/m0*((H01*H31-H03*H11)*(P2*S3-P3*S2)                    &
+ (H01*H32-H03*H12)*(P3*S1-P1*S3) + (H01*H33-H03*H13)*(P1*S2-P2*S1)   & 
+ (H11*H32-H13*H12)*(P0*S3-P3*S0) + (H11*H33-H13*H13)*(P2*S0-P0*S2)   &
+ (H21*H33-H23*H13)*(P0*S1-P1*S0))
S_23 = eta0123/m0*((H02*H31-H03*H21)*(P2*S3-P3*S2)                    &
+ (H02*H32-H03*H22)*(P3*S1-P1*S3) + (H02*H33-H03*H23)*(P1*S2-P2*S1)   & 
+ (H12*H32-H13*H22)*(P0*S3-P3*S0) + (H12*H33-H13*H23)*(P2*S0-P0*S2)   &
+ (H22*H33-H23*H23)*(P0*S1-P1*S0))




end subroutine calculate_spintensor

subroutine calculate_FourVelocity(PVector, metric, Xprime)

!Arguments
real(kind=dp), dimension(4), intent(in) :: PVector
real(kind=dp), dimension(4,4), intent(in) :: metric
real(kind=dp), dimension(4), intent(out) :: Xprime

!Other
real(kind=dp) :: M, delta, RPS_0, RPS_1, RPS_2, RPS_3
real(kind=dp) :: P0, P1, P2, P3
real(kind=dp) :: Vsq, PV
real(kind=dp) :: g00,g11,g22,g33,g30


g00 = metric(1,1)
g11 = metric(2,2)
g22 = metric(3,3)
g33 = metric(4,4)
g30 = metric(4,1)





!read in relevant data
P0 = PVector(1)
P1 = PVector(2)
P2 = PVector(3)
P3 = PVector(4)


M = 1.0_dp !legacy
 
delta = 2.0_dp*(m0**2 + lambda*                                                             &
(R_0101*S_01**2+R_0202*S_02**2+R_2323*S_23**2+R_0203*S_02*S_03+R_1301*S_13*S_01             & 
+R_1201*S_12*S_01+R_0313*S_03*S_13+R_0312*S_03*S_12+R_0113*S_01*S_13+R_0112*S_01*S_12       &
+R_0301*S_03*S_01+R_1212*S_12**2+R_0201*S_02*S_01+R_0303*S_03**2+R_0302*S_03*S_02           &
+R_1312*S_13*S_12+R_0103*S_01*S_03+R_0102*S_01*S_02+R_1313*S_13**2+R_0123*S_01*S_23         &
+R_1202*S_12*S_02+R_1203*S_12*S_03+R_2313*S_23*S_13+R_2312*S_23*S_12+R_1302*S_13*S_02       &
+R_2303*S_23*S_03+R_2302*S_23*S_02+R_0213*S_02*S_13+R_0212*S_02*S_12+R_2301*S_23*S_01       &
+R_1303*S_13*S_03+R_1213*S_12*S_13+R_1223*S_12*S_23+R_0223*S_02*S_23+R_1323*S_13*S_23       &
+R_0323*S_03*S_23)/4.0_dp)


RPS_0 = S_01**2*R_1301*P3+S_01**2*R_1201*P2+S_02*R_2312*P3*S_12+S_01*R_1203*P2*S_03         &
+S_02*R_2313*P3*S_13-S_01**2*R_0101*P0+S_01*R_1213*P2*S_13-S_02*R_0203*P0*S_03-S_02*R_1201  &
 *P1*S_01-S_03*R_2301*P2*S_01-S_03*R_1301*P1*S_01-S_03*R_0301*P0*S_01-S_02*R_0201*P0*S_01    &
-S_02*R_1203*P1*S_03-S_01*R_0103*P0*S_03-S_03*R_2302*P2*S_02-S_03*R_1302*P1*S_02-S_01       &
 *R_0112*P0*S_12-S_03*R_2312*P2*S_12-S_03*R_1312*P1*S_12-S_03*R_0312*P0*S_12-S_02*R_1212*P1  &
 *S_12-S_02*R_0212*P0*S_12-S_01*R_0102*P0*S_02-S_01*R_0113*P0*S_13-S_03*R_2313*P2*S_13       &
-S_03*R_1313*P1*S_13-S_03*R_0313*P0*S_13-S_02*R_1213*P1*S_13-S_02*R_0213*P0*S_13-S_02       &
 *R_1223*P1*S_23-S_01*R_0123*P0*S_23-S_03*R_2323*P2*S_23-S_03*R_1323*P1*S_23+S_01*R_1302*P3  &
 *S_02+S_01*R_1202*P2*S_02+S_01*R_1313*P3*S_13+S_01*R_1312*P3*S_12+S_01*R_1212*P2*S_12       &
+S_02**2*R_2302*P3-S_03*R_0302*P0*S_02-S_02**2*R_0202*P0-S_02**2*R_1202*P1-S_03**2*R_0303   &
 *P0-S_03**2*R_1303*P1-S_03**2*R_2303*P2+S_02*R_2323*P3*S_23-S_03*R_0323*P0*S_23+S_01*R_1323 &
 *P3*S_23+S_01*R_1223*P2*S_23-S_02*R_0223*P0*S_23+S_02*R_2301*P3*S_01+S_01*R_1303*P3*S_03    &
+S_02*R_2303*P3*S_03

RPS_1 = S_12**2*R_2312*P3+S_12*R_2303*P3*S_03+S_12*R_2313*P3*S_13+S_12*R_2323*P3*S_23       &
-S_01**2*R_0101*P1-S_01**2*R_0201*P2-S_01**2*R_0301*P3-S_12**2*R_0212*P0-S_13**2*R_1313*P1  &
-S_13**2*R_0313*P0-S_01*R_0103*P1*S_03-S_13*R_2301*P2*S_01-S_13*R_0301*P0*S_01-S_13*R_1312  &
 *P1*S_12-S_13*R_0312*P0*S_12-S_13*R_2312*P2*S_12-S_12*R_0201*P0*S_01-S_13*R_2302*P2*S_02    &
-S_13*R_1302*P1*S_02-S_12*R_1203*P1*S_03-S_12*R_0203*P0*S_03-S_12*R_1201*P1*S_01-S_01       &
 *R_0303*P3*S_03-S_01*R_0203*P2*S_03-S_12**2*R_1212*P1-S_01*R_0323*P3*S_23-S_01*R_0223*P2    &
 *S_23-S_01*R_0123*P1*S_23-S_12*R_1223*P1*S_23-S_12*R_0223*P0*S_23-S_13*R_0323*P0*S_23       &
-S_13*R_2323*P2*S_23-S_13*R_1323*P1*S_23+S_12*R_2302*P3*S_02-S_12*R_1202*P1*S_02-S_01       &
 *R_0202*P2*S_02-S_13*R_0302*P0*S_02-S_01*R_0113*P1*S_13-S_01*R_0302*P3*S_02-S_12*R_0202     &
 *P0*S_02-S_01*R_0102*P1*S_02-S_12*R_0213*P0*S_13-S_01*R_0313*P3*S_13-S_01*R_0213*P2*S_13    &
-S_13**2*R_2313*P2-S_01*R_0212*P2*S_12-S_01*R_0112*P1*S_12-S_12*R_1213*P1*S_13-S_01*R_0312  &
 *P3*S_12-S_13*R_1303*P1*S_03-S_13*R_0303*P0*S_03-S_13*R_1301*P1*S_01-S_13*R_2303*P2*S_03    &
+S_12*R_2301*P3*S_01

RPS_2 = -S_12*R_1201*P2*S_01-S_02*R_0301*P3*S_01-S_02*R_0201*P2*S_01-S_23*R_1312*P1*S_12    &
-S_02*R_0101*P1*S_01-S_02*R_0112*P1*S_12-S_02*R_0312*P3*S_12+S_12**2*R_0112*P0-S_12*R_1203  &
 *P2*S_03+S_12*R_0101*P0*S_01-S_23*R_0312*P0*S_12+S_12*R_0103*P0*S_03-S_23*R_2312*P2*S_12    &
+S_12*R_0113*P0*S_13-S_02**2*R_0102*P1-S_12*R_1302*P3*S_02-S_23*R_1303*P1*S_03+S_12*R_0102  &
 *P0*S_02-S_02**2*R_0302*P3-S_23*R_0303*P0*S_03-S_12*R_1313*P3*S_13-S_12**2*R_1312*P3        &
-S_12**2*R_1212*P2-S_23**2*R_0323*P0-S_23*R_2303*P2*S_03-S_23**2*R_1323*P1-S_02**2*R_0202   &
 *P2-S_23*R_0301*P0*S_01-S_12*R_1301*P3*S_01-S_23*R_2301*P2*S_01-S_23*R_1301*P1*S_01-S_02    &
 *R_0212*P2*S_12-S_12*R_1303*P3*S_03-S_02*R_0323*P3*S_23-S_23**2*R_2323*P2-S_02*R_0303*P3    &
 *S_03+S_12*R_0123*P0*S_23-S_12*R_1223*P2*S_23-S_12*R_1213*P2*S_13-S_02*R_0203*P2*S_03-S_23  &
 *R_2302*P2*S_02-S_02*R_0103*P1*S_03-S_23*R_1302*P1*S_02-S_23*R_0302*P0*S_02-S_23*R_2313*P2  &
 *S_13-S_23*R_1313*P1*S_13-S_02*R_0113*P1*S_13-S_02*R_0313*P3*S_13-S_02*R_0213*P2*S_13       &
-S_12*R_1202*P2*S_02-S_12*R_1323*P3*S_23-S_23*R_0313*P0*S_13-S_02*R_0223*P2*S_23-S_02       &
 *R_0123*P1*S_23

RPS_3 = S_23*R_0201*P0*S_01-S_13*R_1301*P3*S_01-S_03*R_0201*P2*S_01-S_03*R_0101*P1*S_01     &
-S_13*R_1201*P2*S_01-S_03*R_0301*P3*S_01-S_23*R_2303*P3*S_03-S_13*R_1303*P3*S_03+S_23       &
 *R_1212*P1*S_12+S_13*R_0103*P0*S_03+S_23*R_0212*P0*S_12-S_23*R_2312*P3*S_12-S_13*R_1203*P2  &
 *S_03+S_13*R_0123*P0*S_23-S_13*R_1323*P3*S_23-S_13*R_1223*P2*S_23-S_03*R_0323*P3*S_23       &
-S_03*R_0123*P1*S_23-S_03*R_0223*P2*S_23+S_23**2*R_0223*P0+S_23**2*R_1223*P1-S_23*R_2302*P3 &
 *S_02+S_13*R_0102*P0*S_02+S_23*R_1202*P1*S_02-S_23**2*R_2323*P3+S_23*R_0202*P0*S_02-S_13    &
 *R_1202*P2*S_02-S_13*R_1302*P3*S_02-S_03*R_0313*P3*S_13+S_13**2*R_0113*P0-S_03*R_0213*P2    &
 *S_13-S_03*R_0102*P1*S_02+S_23*R_0213*P0*S_13-S_13*R_1312*P3*S_12-S_13*R_1212*P2*S_12-S_03  &
 *R_0302*P3*S_02-S_03*R_0202*P2*S_02-S_03*R_0113*P1*S_13+S_23*R_1213*P1*S_13-S_13**2*R_1213  &
 *P2-S_03*R_0212*P2*S_12-S_03*R_0112*P1*S_12-S_23*R_2313*P3*S_13+S_13*R_0112*P0*S_12-S_13**2 &
 *R_1313*P3-S_03*R_0312*P3*S_12+S_23*R_1203*P1*S_03+S_23*R_0203*P0*S_03+S_13*R_0101*P0*S_01  &
+S_23*R_1201*P1*S_01-S_03**2*R_0303*P3-S_23*R_2301*P3*S_01-S_03**2*R_0203*P2-S_03**2*R_0103 &
 *P1



Xprime(1) = - (P0 + lambda*RPS_0/Delta)/m0**2
Xprime(2) = - (P1 + lambda*RPS_1/Delta)/m0**2
Xprime(3) = - (P2 + lambda*RPS_2/Delta)/m0**2
Xprime(4) = - (P3 + lambda*RPS_3/Delta)/m0**2





Vsq = g00*Xprime(1)**2 + g11*Xprime(2)**2 + g22*Xprime(3)**2 + g33*Xprime(4)**2 + 2.0_dp*g30*Xprime(1)*Xprime(4)





PV = -sqrt(-1.0_dp/Vsq)


Xprime = Xprime * PV



end subroutine calculate_FourVelocity




subroutine calculate_FourMom(VVector, PVector, Pprime)
!Arguments
real(kind=dp), intent(in), dimension(4) :: VVector, PVector !Velocity and momentum vector

real(kind=dp), intent(out), dimension(4) :: Pprime 

!Other
real (kind=dp) :: V0,V1,V2,V3
real(kind=dp) :: P0,P1,P2,P3 
real(kind=dp) :: Pdot1_0, Pdot1_1, Pdot1_2, Pdot1_3
real(kind=dp) :: RS0_0, RS0_1, RS0_2, RS0_3, RS1_0, RS1_1, RS1_2, RS1_3,                     &
                 RS2_0, RS2_1, RS2_2, RS2_3, RS3_0, RS3_1, RS3_2, RS3_3
real(kind=dp) :: Pdot2_0, Pdot2_1, Pdot2_2, Pdot2_3 

V0 = VVector(1)
V1 = VVector(2)
V2 = VVector(3)
V3 = VVector(4)


P0 = PVector(1)
P1 = PVector(2)
P2 = PVector(3)
P3 = PVector(4)


Pdot1_0 = -G0_00*V0*P0-G0_01*V0*P1-G0_02*V0*P2-G0_03*V0*P3-G0_01*V1*P0-G0_11*V1*P1-G0_12*V1 &
 *P2-G0_13*V1*P3-G0_02*V2*P0-G0_12*V2*P1-G0_22*V2*P2-G0_23*V2*P3-G0_03*V3*P0-G0_13*V3*P1     &
-G0_23*V3*P2-G0_33*V3*P3

Pdot1_1 = -G1_00*V0*P0-G1_01*V0*P1-G1_02*V0*P2-G1_03*V0*P3-G1_01*V1*P0-G1_11*V1*P1-G1_12*V1 &
 *P2-G1_13*V1*P3-G1_02*V2*P0-G1_12*V2*P1-G1_22*V2*P2-G1_23*V2*P3-G1_03*V3*P0-G1_13*V3*P1     &
-G1_23*V3*P2-G1_33*V3*P3

Pdot1_2 = -G2_00*V0*P0-G2_01*V0*P1-G2_02*V0*P2-G2_03*V0*P3-G2_01*V1*P0-G2_11*V1*P1-G2_12*V1 &
 *P2-G2_13*V1*P3-G2_02*V2*P0-G2_12*V2*P1-G2_22*V2*P2-G2_23*V2*P3-G2_03*V3*P0-G2_13*V3*P1     &
-G2_23*V3*P2-G2_33*V3*P3

Pdot1_3 = -G3_00*V0*P0-G3_01*V0*P1-G3_02*V0*P2-G3_03*V0*P3-G3_01*V1*P0-G3_11*V1*P1-G3_12*V1 &
 *P2-G3_13*V1*P3-G3_02*V2*P0-G3_12*V2*P1-G3_22*V2*P2-G3_23*V2*P3-G3_03*V3*P0-G3_13*V3*P1     &
-G3_23*V3*P2-G3_33*V3*P3


RS0_0 = R0_001*S_01 + R0_002*S_02 + R0_003*S_03 + R0_012*S_12 + R0_013*S_13 + R0_023*S_23
RS0_1 = R0_101*S_01 + R0_102*S_02 + R0_103*S_03 + R0_112*S_12 + R0_113*S_13 + R0_123*S_23
RS0_2 = R0_201*S_01 + R0_202*S_02 + R0_203*S_03 + R0_212*S_12 + R0_213*S_13 + R0_223*S_23
RS0_3 = R0_301*S_01 + R0_302*S_02 + R0_303*S_03 + R0_312*S_12 + R0_313*S_13 + R0_323*S_23

RS1_0 = R1_001*S_01 + R1_002*S_02 + R1_003*S_03 + R1_012*S_12 + R1_013*S_13 + R1_023*S_23
RS1_1 = R1_101*S_01 + R1_102*S_02 + R1_103*S_03 + R1_112*S_12 + R1_113*S_13 + R1_123*S_23
RS1_2 = R1_201*S_01 + R1_202*S_02 + R1_203*S_03 + R1_212*S_12 + R1_213*S_13 + R1_223*S_23
RS1_3 = R1_301*S_01 + R1_302*S_02 + R1_303*S_03 + R1_312*S_12 + R1_313*S_13 + R1_323*S_23

RS2_0 = R2_001*S_01 + R2_002*S_02 + R2_003*S_03 + R2_012*S_12 + R2_013*S_13 + R2_023*S_23
RS2_1 = R2_101*S_01 + R2_102*S_02 + R2_103*S_03 + R2_112*S_12 + R2_113*S_13 + R2_123*S_23
RS2_2 = R2_201*S_01 + R2_202*S_02 + R2_203*S_03 + R2_212*S_12 + R2_213*S_13 + R2_223*S_23
RS2_3 = R2_301*S_01 + R2_302*S_02 + R2_303*S_03 + R2_312*S_12 + R2_313*S_13 + R2_323*S_23

RS3_0 = R3_001*S_01 + R3_002*S_02 + R3_003*S_03 + R3_012*S_12 + R3_013*S_13 + R3_023*S_23
RS3_1 = R3_101*S_01 + R3_102*S_02 + R3_103*S_03 + R3_112*S_12 + R3_113*S_13 + R3_123*S_23
RS3_2 = R3_201*S_01 + R3_202*S_02 + R3_203*S_03 + R3_212*S_12 + R3_213*S_13 + R3_223*S_23
RS3_3 = R3_301*S_01 + R3_302*S_02 + R3_303*S_03 + R3_312*S_12 + R3_313*S_13 + R3_323*S_23



Pdot2_0 = -(RS0_0*V0 + RS0_1*V1 + RS0_2*V2 + RS0_3*V3)

Pdot2_1 = -(RS1_0*V0 + RS1_1*V1 + RS1_2*V2 + RS1_3*V3)

Pdot2_2 = -(RS2_0*V0 + RS2_1*V1 + RS2_2*V2 + RS2_3*V3)

Pdot2_3 = -(RS3_0*V0 + RS3_1*V1 + RS3_2*V2 + RS3_3*V3)


Pprime(1) = Pdot1_0 + lambda*Pdot2_0
Pprime(2) = Pdot1_1 + lambda*Pdot2_1
Pprime(3) = Pdot1_2 + lambda*Pdot2_2
Pprime(4) = Pdot1_3 + lambda*Pdot2_3






end subroutine calculate_FourMom


subroutine calculate_FourSpin(VVector,PVector,SVector, Sprime)
!Arguments
real(kind=dp), dimension(4), intent(in) :: VVector, PVector,SVector
real(kind=dp), dimension(4), intent(out) :: Sprime
!Other
real (kind=dp) :: V0,V1,V2,V3
real(kind=dp) :: S0,S1,S2,S3 
real(kind=dp) :: P0,P1,P2,P3
real(kind=dp) :: Scoeff
real(kind=dp) :: Sdot2_0, Sdot2_1, Sdot2_2, Sdot2_3
real(kind=dp) :: Sdot1_0, Sdot1_1, Sdot1_2, Sdot1_3



V0 = VVector(1)
V1 = VVector(2)
V2 = VVector(3)
V3 = VVector(4)


P0 = PVector(1)
P1 = PVector(2)
P2 = PVector(3)
P3 = PVector(4)


S0 = SVector(1)
S1 = SVector(2)
S2 = SVector(3)
S3 = SVector(4)



Scoeff = -(R_0212*S0*V2*S_12-R_1213*S2*V1*S_13-R_1212*S2*V1*S_12-R_1201*S2*V1*S_01-R_2302   &
 *S3*V2*S_02+R_0312*S0*V3*S_12+R_0313*S0*V3*S_13+R_0302*S0*V3*S_02-R_0323*S3*V0*S_23+R_2301  &
 *S2*V3*S_01-R_1303*S3*V1*S_03-R_1202*S2*V1*S_02-R_1203*S2*V1*S_03-R_1223*S2*V1*S_23-R_0301  &
 *S3*V0*S_01-R_0312*S3*V0*S_12-R_1302*S3*V1*S_02-R_1301*S3*V1*S_01+R_2323*S2*V3*S_23-R_0303  &
 *S3*V0*S_03-R_2323*S3*V2*S_23-R_0313*S3*V0*S_13-R_0302*S3*V0*S_02-R_2312*S3*V2*S_12-R_2301  &
 *S3*V2*S_01-R_2313*S3*V2*S_13+R_0113*S0*V1*S_13+R_0213*S0*V2*S_13+R_0223*S0*V2*S_23+R_0112  &
 *S0*V1*S_12-R_0102*S1*V0*S_02+R_0203*S0*V2*S_03-R_0112*S1*V0*S_12-R_0123*S1*V0*S_23+R_0103  &
 *S0*V1*S_03+R_0102*S0*V1*S_02-R_1313*S3*V1*S_13-R_1312*S3*V1*S_12-R_1323*S3*V1*S_23-R_2303  &
 *S3*V2*S_03-R_0113*S1*V0*S_13-R_0101*S1*V0*S_01-R_0103*S1*V0*S_03+R_1223*S1*V2*S_23+R_1213  &
 *S1*V2*S_13+R_1302*S1*V3*S_02+R_0303*S0*V3*S_03+R_0301*S0*V3*S_01+R_0323*S0*V3*S_23-R_0223  &
 *S2*V0*S_23+R_1203*S1*V2*S_03+R_1202*S1*V2*S_02+R_1201*S1*V2*S_01+R_1323*S1*V3*S_23+R_1303  &
 *S1*V3*S_03+R_1312*S1*V3*S_12+R_1212*S1*V2*S_12-R_0213*S2*V0*S_13-R_0201*S2*V0*S_01-R_0202  &
 *S2*V0*S_02-R_0212*S2*V0*S_12+R_1313*S1*V3*S_13+R_1301*S1*V3*S_01-R_0203*S2*V0*S_03+R_2303  &
 *S2*V3*S_03+R_0101*S0*V1*S_01+R_2302*S2*V3*S_02+R_2313*S2*V3*S_13+R_0123*S0*V1*S_23+R_2312  &
 *S2*V3*S_12+R_0202*S0*V2*S_02+R_0201*S0*V2*S_01)/m0**2



Sdot2_0 = Scoeff*P0
Sdot2_1 = Scoeff*P1
Sdot2_2 = Scoeff*P2
Sdot2_3 = Scoeff*P3


Sdot1_0 = -G0_00*V0*S0-G0_01*V0*S1-G0_02*V0*S2-G0_03*V0*S3-G0_01*V1*S0-G0_11*V1*S1-G0_12*V1 &
 *S2-G0_13*V1*S3-G0_02*V2*S0-G0_12*V2*S1-G0_22*V2*S2-G0_23*V2*S3-G0_03*V3*S0-G0_13*V3*S1     &
-G0_23*V3*S2-G0_33*V3*S3

Sdot1_1 = -G1_00*V0*S0-G1_01*V0*S1-G1_02*V0*S2-G1_03*V0*S3-G1_01*V1*S0-G1_11*V1*S1-G1_12*V1 &
 *S2-G1_13*V1*S3-G1_02*V2*S0-G1_12*V2*S1-G1_22*V2*S2-G1_23*V2*S3-G1_03*V3*S0-G1_13*V3*S1     &
-G1_23*V3*S2-G1_33*V3*S3

Sdot1_2 = -G2_00*V0*S0-G2_01*V0*S1-G2_02*V0*S2-G2_03*V0*S3-G2_01*V1*S0-G2_11*V1*S1-G2_12*V1 &
 *S2-G2_13*V1*S3-G2_02*V2*S0-G2_12*V2*S1-G2_22*V2*S2-G2_23*V2*S3-G2_03*V3*S0-G2_13*V3*S1     &
-G2_23*V3*S2-G2_33*V3*S3

Sdot1_3 = -G3_00*V0*S0-G3_01*V0*S1-G3_02*V0*S2-G3_03*V0*S3-G3_01*V1*S0-G3_11*V1*S1-G3_12*V1 &
 *S2-G3_13*V1*S3-G3_02*V2*S0-G3_12*V2*S1-G3_22*V2*S2-G3_23*V2*S3-G3_03*V3*S0-G3_13*V3*S1     &
-G3_23*V3*S2-G3_33*V3*S3



Sprime(1) = Sdot1_0 + lambda*Sdot2_0
Sprime(2) = Sdot1_1 + lambda*Sdot2_1
Sprime(3) = Sdot1_2 + lambda*Sdot2_2
Sprime(4) = Sdot1_3 + lambda*Sdot2_3


end subroutine calculate_FourSpin


!Properties of the plasma
subroutine plasma_fr(r,fr)
real(kind=dp), intent(in) :: r
real(kind=dp), intent(out) :: fr
fr = r**0.90_dp 
end subroutine plasma_fr


subroutine plasma_fr_deriv(r,fr)
real(kind=dp), intent(in) :: r
real(kind=dp), intent(out) :: fr
fr = 0.90_dp*r**(-0.10_dp)
end subroutine plasma_fr_deriv



subroutine plasma_ft(r,ft)
real(kind=dp), intent(in) :: r
real(kind=dp), intent(out) :: ft
ft = 0.0_dp
end subroutine plasma_ft


subroutine plasma_ft_deriv(r,ft)
real(kind=dp), intent(in) :: r
real(kind=dp), intent(out) :: ft
ft = 0.0_dp
end subroutine plasma_ft_deriv






end module metric
