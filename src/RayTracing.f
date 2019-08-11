module RayTracing


use parameters
use constants
use metric
use NumericalMethods
use derivatives

implicit none

public RT

private initial_conditions

contains

subroutine RT()
integer(kind=dp) :: stat,i

real(kind=dp), dimension(:,:), allocatable :: MPDData
real(kind=dp) :: t,r,theta,phi,tau
real(kind=dp), dimension(4) :: xvector, uvector, svector
real(kind=dp) :: p1,p2,p3,p4
real(kind=dp) :: mag, rstart
real(kind=dp), dimension(6) :: v
real(kind=dp), dimension(4) :: c
real(kind=dp),dimension(:,:),allocatable :: PlotData !Big array to save all data
integer(kind=dp) :: counter , jj
real(kind=dp) :: xC, yC, zC,cpt
real(kind=dp),dimension(4,4) :: metric_covar
real(kind=dp) :: gtt, gpp, gtp, aa,bb,cc, Eobs
character(len=200) :: FileName,IDStr
integer(kind=dp) :: FileID



allocate(MPDData(nrows,entries+1))
allocate(PlotData(int(1e5),6))



!Load the MPD Data file
open(unit=10, file=MPDBinaryData , form='unformatted',access='stream')

!The file is composed of multiple arrays of size nrows x entries+1.
!We load each array in turn.

stat = 0
do while (stat .EQ. 0)


    !For each array
    read(10,iostat=stat) MPDData

    !For each row in array   Bring OpenMp in here
    do i = 1,nrows

    if (MPDData(i,2) .EQ. 0.0_dp) then
    exit
    endif




    !Get 4 velocity
    uvector(1) = MPDData(i,5)/m0
    uvector(2) = MPDData(i,6)/m0
    uvector(3) = MPDData(i,7)/m0
    uvector(4) = MPDData(i,8)/m0

    !Get position coords
    xvector(1) = MPDData(i,13) !tau
    xvector(2:4) = MPDData(i,2:4)

    !Get spin
    svector = MPDData(i,9:12)
    
 
    
    !Check that the 4-velocity satisfies mag = -1
 
    r = MPDData(i,2)
    theta = MPDData(i,3)
    call mag_4(r,theta,uvector,mag)
    !print *, 'Magnitude of 4-velocity = ', mag
    
    
    !Set the energy at +infinity
    Eobs = 2.0_dp*PI*1.0d9

    !SEt initial conditions
    !Vector v = (r,that,phi,t,pr,ptheta) c = constants + stepsize
    call initial_conditions(xvector,uvector,svector,Eobs, &
                            v,c)
    
    rstart = v(1)
    counter = 1
    do while (v(1) .GT. Rhor .and. v(1) .LT. rstart*10.0_dp)
    call rk_geodesic(v,c)
    
 !   print *, counter, v(1)


    if (plot_RT .EQ. 1) then
    PlotData(counter,:) = v
    endif


    counter = counter + 1
    enddo



    !Save to file

    
FileID = i*100
write (IDStr, "(I5)") FileID
FileName = trim(adjustl(RTPath))//'RT_'//trim(adjustl(IDStr))//'.txt'

open(unit = FileID, file = FileName, status = 'replace', form='formatted')



do jj=1,counter-1
    xC = sqrt(PlotData(jj,1)**2 + a**2) * sin(PlotData(jj,2))*cos(PlotData(jj,3))
    yC = sqrt(PlotData(jj,1)**2 + a**2) * sin(PlotData(jj,2))*sin(PlotData(jj,3))
    zC = PlotData(jj,1) * cos(PlotData(jj,2))
    
 !   print *, xC, yC, zC

    write(FileID,*) xC,yC,zC

enddo


close(FileID)





    enddo

enddo









close(10)

end subroutine RT






subroutine initial_conditions(xi,ui,si,Eobs, &
                              v,c)
!Arguments

real(kind=dp), dimension(4), intent(in) :: xi !position
real(kind=dp), dimension(4), intent(in) :: ui !4 velocity
real(kind=dp), dimension(4), intent(in) :: si !4 spin
real(kind=dp),intent(in) :: Eobs !The observation 'Frequency'
real(kind=dp), dimension(6), intent(out) :: v
real(kind=dp), dimension(4), intent(out) :: c

!Other
real(kind=dp),dimension(4) :: ki
real(kind=dp) :: psi,tau, s1,s2,s3,r,theta,phi
real(kind=dp) :: sx,sy,sz, stheta, sphi
real(kind=dp),dimension(3) :: ki_rot, ki_space
real(kind=dp), dimension(3,3) :: Rz, Ry
real(kind=dp) :: r_dot, theta_dot, phi_dot,delta
real(kind=dp) :: fr,ft, omega2, sigma, Enorm,Eprime, Eprime2, E2
real(kind=dp) :: pr, ptheta, Lz, B2, kappa, mm
real(kind=dp) :: xdot, ydot,zdot
real(kind=dp) :: st, sp !stheta, sphi
!NOTE: CAN WE OPETATE ON THE WHOLE ARRAY TO SPEED UP? 

r = xi(2)
theta = xi(3)
phi = xi(4)

sigma = r**2 +a**2*cos(theta)**2
delta = r**2 -2.0_dp*r + a**2

!Set dipole axis
tau = xi(1)
psi = 2.0_dp*PI*tau/SpinPeriod

!Define direction tangent vector
ki(1) = 1.0_dp
ki(2) = sin(chi)*cos(psi)
ki(3) = sin(chi)*sin(psi)
ki(4) = cos(chi)


!Rotate it to account for precession of spin axis
s1 = si(2)
s2 = si(3)
s3 = si(4)

!no need for BL here - defined locally. I think
sx = s1*sin(theta)*cos(phi) + s2*r*cos(theta)*cos(phi) - s3*r*sin(theta)*sin(phi)
sy = s1*sin(theta)*sin(phi) + s2*r*cos(theta)*sin(phi) + s3*r*sin(theta)*cos(phi)
sz = s1*cos(theta) - s2*r*sin(theta)

st = atan2(sqrt(sx**2 + sy**2),sz)
sp = atan2(sy,sx)



Rz = 0.0_dp
Ry= 0.0_dp
Rz(1,1) = cos(sp)
Rz(1,2) = -sin(sp)
Rz(2,1) = sin(sp)
Rz(2,2) = cos(sp)
Rz(3,3) = 1.0_dp

Ry(1,1) = cos(st)
Ry(1,3) = sin(st)
Ry(3,1) = -sin(st)
Ry(3,3) = cos(st)
Ry(2,2) = 1.0_dp





ki_space = ki(2:4)
ki_rot  = MATMUL(Rz,MATMUL(Ry,ki_space))

ki(2:4) = ki_space



mm = sqrt(r**2 + a**2)



!ki is still xdot ,ydoy,zdot 
!need a transform to rdot,zdot,thetadot in tetrad frame?


!This is the cartesian tetrad frame
xdot = ki(2)
ydot = ki(3)
zdot = ki(4)

!This is the polar tetrad frame





ki(2) = sin(theta)*cos(phi) * xdot + &
        sin(theta)*cos(phi) * ydot + &
        cos(theta) * zdot


ki(3) = cos(theta)*cos(phi) * xdot + &
        cos(theta)*sin(phi) * ydot - &
        sin(theta)


ki(4) = -sin(phi)*xdot + cos(phi)*ydot



!mm = r
r_dot = mm*r*sin(theta)*cos(phi)*xdot/sigma &
       +mm*r*sin(theta)*sin(phi)*ydot/sigma &
       +mm**2*cos(theta)*zdot/sigma



theta_dot = (mm*cos(theta)*cos(phi) * xdot &
           +mm*cos(theta)*sin(phi) * ydot &
           -r*sin(theta)* zdot&
           )/sigma


phi_dot = (-sin(phi)*xdot + cos(phi)*ydot)/(mm*sin(theta))


call transform_to_global(xi,ui,ki)
r_dot = ki(2)
theta_dot = ki(3)
phi_dot = ki(4)



!Define plasma freq

B2 =  N*4.0_dp*PI*electron_charge**2 / electron_mass

call plasma_fr(r,fr)
call plasma_ft(r,ft)
omega2 = B2 * (fr+ft)/sigma


!Define E normalizatiopn
Enorm = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2) + delta*(sin(theta)*phi_dot)**2
Eprime2 = (Eobs**2 - (sigma-2.0_dp*r)*omega2/sigma)/Enorm
Eprime = sqrt(Eprime2)


!Correct for E normalisation
r_dot = r_dot * Eprime
theta_dot = theta_dot*Eprime
phi_dot = phi_dot*Eprime



!Check energies are equal
E2 = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2 + omega2/sigma) + delta*(sin(theta)*phi_dot)**2

!print *, 'Energy check:', Eobs/ sqrt(E2)

pr = r_dot * sigma/delta
ptheta = sigma*theta_dot
Lz = (sigma*delta*phi_dot - 2.0_dp*a*r*Eobs)*sin(theta)**2 / (sigma-2.0_dp*r)



!Normalise to E = 1
pr = pr/Eobs
ptheta = ptheta/Eobs
Lz = Lz/Eobs
B2 = B2 / Eobs**2

kappa = ptheta**2 + Lz**2/sin(theta)**2 + a**2*sin(theta)**2

v(1) = r
v(2) = theta
v(3) = phi
v(4) = 0.0_dp !t
v(5) = pr
v(6) = ptheta


!Dont declare these globally as need to be careful when running in parallel
c(1) = Lz
c(2) = kappa
c(3) = B2
c(4) = 1.0d-6 !Initial stepsize for RT

end subroutine initial_conditions

subroutine transform_to_global(xi,ui,ki)
!Arguments
real(kind=dp), dimension(4), intent(in) :: xi,ui
real(kind=dp), dimension(4), intent(inout) :: ki

!Other
real(kind=dp) :: r, theta,delta,N1,N2,N3,mag
real(kind=dp), dimension(4,4) :: metric_contra, metric_covar, transform_matrix
real(kind=dp), dimension(4) :: u_covar, u_contra, ki_global
integer(kind=dp) :: i

!Load position
r = xi(2)
theta = xi(3)


!Compute metric
call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)



!Get contra and covar 4-velocity
u_contra = ui
do i = 1,4
u_covar(i) = metric_covar(i,1) * u_contra(1) + &
             metric_covar(i,2) * u_contra(2) + &
             metric_covar(i,3) * u_contra(3) + &
             metric_covar(i,4) * u_contra(4) 
enddo


delta = r**2.0_dp + a**2.0_dp - 2.0_dp*r
N1 = sqrt(- metric_covar(2,2) * (u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4)) * (1.0_dp + u_covar(3)*u_contra(3)) )
N2 = sqrt(metric_covar(3,3) * (1.0_dp + u_covar(3) * u_contra(3)) )
N3 = sqrt(-(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))*delta*sin(theta)**2)



transform_matrix(1,:) = u_contra

transform_matrix(2,1) = u_covar(2)*u_contra(1)/N1 
transform_matrix(2,2) = -(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))/N1
transform_matrix(2,3) = 0.0_dp
transform_matrix(2,4) = u_covar(2)*u_contra(4)/N1


transform_matrix(3,1) = u_covar(3)*u_contra(1)/N2
transform_matrix(3,2) = u_covar(3)*u_contra(2) / N2
transform_matrix(3,3) = (1.0_dp + u_covar(3)*u_contra(3))/N2
transform_matrix(3,4) = u_covar(3)*u_contra(4)/N2



transform_matrix(4,1) = u_covar(4)/N3
transform_matrix(4,2) = 0.0_dp
transform_matrix(4,3) = 0.0_dp
transform_matrix(4,4) = -u_covar(1)/N3



!Note this is not just a MatMul. See Kulkarni et al.



do i = 1,4

ki_global(i) = transform_matrix(1,i)*ki(1) + &
               transform_matrix(2,i)*ki(2) + &
               transform_matrix(3,i)*ki(3) + &
               transform_matrix(4,i)*ki(4) 
enddo


!A check
!call mag_4(r,theta,ki_global,mag)
!print *, 'Magnitude checks'
!print *, -(ki(1))**2 + ki(2)**2 + ki(3)**2 + ki(4)**2
!print *, mag
ki = ki_global

end subroutine transform_to_global

end module RayTracing
