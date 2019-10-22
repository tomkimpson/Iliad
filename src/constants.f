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
real(kind=dp), parameter :: Rhor = 1.0_dp+sqrt(1.0_dp-a**2) + 1.0d-2 !Horizon + eps
real(kind=dp), PARAMETER :: electron_charge = 4.80320425D-10 !CGS
real(kind=dp), PARAMETER :: electron_mass = 9.10938356D-28 !CGS

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
integer(kind=dp), parameter :: entries = 13 !Number of differetnai eqns 4x(position,spin,momentum + proper time)
real(kind=dp), parameter :: hs_natural = hs * convert_s
real(kind=dp), parameter :: SpinPeriod = p0 * convert_s


!MPD I/O
integer(kind=dp), parameter :: nrows = 1d4
character(len=200) :: path, MPDBinaryData, MPDFormatData, RTPath



!Ray tracing



!Numerical integration
!Cash-Karp parameters and integration constants
real(kind = dp) :: B21=1.0_dp/5.0_dp
real(kind = dp) :: B31 = 3.0_dp/40.0_dp , B32 = 9.0_dp/40.0_dp
real(kind = dp) :: B41 = 3.0_dp/10.0_dp, B42 = -9.0_dp/10.0_dp, B43 = 6.0_dp/5.0_dp 
real(kind = dp) :: B51 = -11.0_dp/54.0_dp, B52 = 5.0_dp/2.0_dp, B53 = -70.0_dp/27.0_dp
real(kind = dp) :: B54 = 35.0_dp/27.0_dp
real(kind = dp) :: B61 = 1631.0_dp/55296.0_dp, B62 = 175.0_dp/512.0_dp, B63 = 575.0_dp/13824.0_dp
real(kind = dp) :: B64 = 44275.0_dp/110592.0_dp, B65 = 253.0_dp/4096.0_dp


real(kind = dp) :: c1 = 37.0_dp/378.0_dp, c3 = 250.0_dp/621.0_dp, c4 = 125.0_dp/594.0_dp
real(kind = dp) :: c6=512.0_dp/1771.0_dp
real(kind = dp) :: cbar1 = 2825.0_dp/27648.0_dp, cbar3 = 18575.0_dp/48384.0_dp
real(kind = dp) :: cbar4=13525.0_dp/55296.0_dp, cbar5 = 277.0_dp/14336.0_dp, cbar6 = 1.0_dp/4.0_dp

real(kind = dp) :: escal = 1.0d19 !SHOULD WE CHANGE THIS? SPEED VS ACCURACY....
real(kind = dp) :: PSHRINK = -0.25_dp,PGROW=-0.20_dp, S=0.90_dp
real(kind=dp), parameter :: dx_eps = 1.0e-22_dp !How small do we push this param?


!Some globally defined MPD variables
real(kind=dp) :: m_sq, s_sq ! mass sqaures and s squared from initial condiitons module
real(kind=dp) :: tau !Proper time. Related to the stepsize

!Christoffel sybols
real(kind=dp) ::    G0_00, G0_01, G0_02, G0_03, G0_11, G0_12, G0_13, G0_22, G0_23, G0_33,       &
                    G1_00, G1_01, G1_02, G1_03, G1_11, G1_12, G1_13, G1_22, G1_23, G1_33,       &
                    G2_00, G2_01, G2_02, G2_03, G2_11, G2_12, G2_13, G2_22, G2_23, G2_33,       &
                    G3_00, G3_01, G3_02, G3_03, G3_11, G3_12, G3_13, G3_22, G3_23, G3_33


!Riemann
REAL(KIND=dp) ::                                                                &
    R_0101, R_0102, R_0103, R_0112, R_0113, R_0123,                             &
    R_0201, R_0202, R_0203, R_0212, R_0213, R_0223,                             &
    R_0301, R_0302, R_0303, R_0312, R_0313, R_0323,                             &    
    R_1201, R_1202, R_1203, R_1212, R_1213, R_1223,                             &
    R_1301, R_1302, R_1303, R_1312, R_1313, R_1323,                             &
    R_2301, R_2302, R_2303, R_2312, R_2313, R_2323,                             &
! 
    R0_001, R0_002, R0_003, R0_012, R0_013, R0_023,                             &
    R0_101, R0_102, R0_103, R0_112, R0_113, R0_123,                             &
    R0_201, R0_202, R0_203, R0_212, R0_213, R0_223,                             &
    R0_301, R0_302, R0_303, R0_312, R0_313, R0_323,                             &
!
    R1_001, R1_002, R1_003, R1_012, R1_013, R1_023,                             &
    R1_101, R1_102, R1_103, R1_112, R1_113, R1_123,                             &
    R1_201, R1_202, R1_203, R1_212, R1_213, R1_223,                             &
    R1_301, R1_302, R1_303, R1_312, R1_313, R1_323,                             &
!
    R2_001, R2_002, R2_003, R2_012, R2_013, R2_023,                             &
    R2_101, R2_102, R2_103, R2_112, R2_113, R2_123,                             &
    R2_201, R2_202, R2_203, R2_212, R2_213, R2_223,                             &
    R2_301, R2_302, R2_303, R2_312, R2_313, R2_323,                             &
!
    R3_001, R3_002, R3_003, R3_012, R3_013, R3_023,                             &
    R3_101, R3_102, R3_103, R3_112, R3_113, R3_123,                             &
    R3_201, R3_202, R3_203, R3_212, R3_213, R3_223,                             &
    R3_301, R3_302, R3_303, R3_312, R3_313, R3_323


!Spin tensor

real(kind=dp) :: S_01, S_02, S_03, S_12, S_13, S_23




end module constants
