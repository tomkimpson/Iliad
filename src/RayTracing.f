module RayTracing


use parameters
use constants
use metric
use NumericalMethods
use derivatives
use IO

implicit none

public RT

private initial_conditions, RT_Forward, RT_Backward, GeneralInitialConditions, shoot, aim

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
    if (stat .NE. 0) then
    exit
    endif

    print *, 'Loaded MPD DATA', stat



    do i = 1,nrows


    if (MPDData(i,2) .EQ. 0.0_dp) then
    !End of data array
    exit
    endif


    if ( mod(real(i,kind=dp), RepRes) .EQ. 0.0_dp) then
 

    !Get the data in terms of vectors
    !Get position coords
    xvector(1) = MPDData(i,13) !tau
    xvector(2:4) = MPDData(i,2:4)









    if (RayTracingDirection .EQ. +1.0_dp) then
    !Forward
    call RT_Forward()
    else if (RayTracingDirection .EQ. -1.0_dp) then
    !Backward
    call RT_Backward(xvector)
    else
    !Exit
    print *, 'Error: You need to set the ray tracing method'
    stop


    endif

            print *, RayTracingDirection
    stop










    endif



    enddo


    stop



















    !For each row in array   Bring OpenMp in here
    do i = 1,nrows


    print *, 'row:', i



    

    !Get 4 velocity
    uvector(1) = MPDData(i,5)/m0
    uvector(2) = MPDData(i,6)/m0
    uvector(3) = MPDData(i,7)/m0
    uvector(4) = MPDData(i,8)/m0

    print *, 'Uvector=', uvector


    !Get spin
    svector = MPDData(i,9:12)
    
 
    
    !Check that the 4-velocity satisfies mag = -1
 
    r = MPDData(i,2)
    theta = MPDData(i,3)
    call mag_4(r,theta,uvector,mag)
    print *, 'Magnitude of 4-velocity = ', mag
    
    
    !Set the energy at +infinity
    Eobs = 2.0_dp*PI*1.0d9

    !SEt initial conditions
    !Vector v = (r,that,phi,t,pr,ptheta) c = constants + stepsize
    
    call initial_conditions(xvector,uvector,svector,Eobs, &
                            v,c)
    
    rstart = v(1)
    counter = 1
 !   do while (v(1) .GT. Rhor .and. v(1) .LT. rstart*10.0_dp)
  
 do while (counter .lt. 10)
! call rk_geodesic(v,c)
    
    print *, counter, v(1),v(2), c(4)


    if (plot_RT .EQ. 1) then
    PlotData(counter,:) = v
    endif


    counter = counter + 1


    enddo


stop

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






!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------

!!Set dipole axis
!tau = xi(1)
!psi = 2.0_dp*PI*tau/SpinPeriod
!
!!Define direction tangent vector
!ki(1) = 1.0_dp
!ki(2) = sin(chi)*cos(psi)
!ki(3) = sin(chi)*sin(psi)
!ki(4) = cos(chi)
!
!
!
!
!!Rotate it to account for precession of spin axis
!s1 = si(2)
!s2 = si(3)
!s3 = si(4)
!
!
!
!!no need for BL here - defined locally. 
!sx = s1*sin(theta)*cos(phi) + s2*r*cos(theta)*cos(phi) - s3*r*sin(theta)*sin(phi)
!sy = s1*sin(theta)*sin(phi) + s2*r*cos(theta)*sin(phi) + s3*r*sin(theta)*cos(phi)
!sz = s1*cos(theta) - s2*r*sin(theta)
!
!st = atan2(sqrt(sx**2 + sy**2),sz)
!sp = atan2(sy,sx)
!
!
!
!Rz = 0.0_dp
!Ry= 0.0_dp
!Rz(1,1) = cos(sp)
!Rz(1,2) = -sin(sp)
!Rz(2,1) = sin(sp)
!Rz(2,2) = cos(sp)
!Rz(3,3) = 1.0_dp
!
!Ry(1,1) = cos(st)
!Ry(1,3) = sin(st)
!Ry(3,1) = -sin(st)
!Ry(3,3) = cos(st)
!Ry(2,2) = 1.0_dp
!
!
!
!
!
!ki_space = ki(2:4)
!ki_rot  = MATMUL(Rz,MATMUL(Ry,ki_space))
!
!ki(2:4) = ki_space
!
!
!
!mm = sqrt(r**2 + a**2)
!
!
!
!!ki is still xdot ,ydoy,zdot 
!!need a transform to rdot,zdot,thetadot in tetrad frame?
!
!
!!This is the cartesian tetrad frame
!xdot = ki(2)
!ydot = ki(3)
!zdot = ki(4)
!
!!This is the polar tetrad frame
!
!
!
!
!
!ki(2) = sin(theta)*cos(phi) * xdot + &
!        sin(theta)*cos(phi) * ydot + &
!        cos(theta) * zdot
!
!
!ki(3) = cos(theta)*cos(phi) * xdot + &
!        cos(theta)*sin(phi) * ydot - &
!        sin(theta)
!
!
!ki(4) = -sin(phi)*xdot + cos(phi)*ydot
!
!
!
!print *, 'ki=', ki
!print *,'xdot', xdot, ydot, zdot
!




!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------
!-----------------------------------------------



!mm = r
r_dot = mm*r*sin(theta)*cos(phi)*xdot/sigma &
       +mm*r*sin(theta)*sin(phi)*ydot/sigma &
       +mm**2*cos(theta)*zdot/sigma



theta_dot = (mm*cos(theta)*cos(phi) * xdot &
           +mm*cos(theta)*sin(phi) * ydot &
           -r*sin(theta)* zdot&
           )/sigma


phi_dot = (-sin(phi)*xdot + cos(phi)*ydot)/(mm*sin(theta))



print *, 'rdot', r_dot, theta_dot, phi_dot
call transform_to_global(xi,ui,ki)
r_dot = ki(2)
theta_dot = ki(3)
phi_dot = ki(4)


print *, 'rdot2', r_dot, theta_dot, phi_dot




!Define plasma freq

B2 =  N*4.0_dp*PI*electron_charge**2 / electron_mass

call plasma_fr(r,fr)
call plasma_ft(r,ft)
omega2 = B2 * (fr+ft)/sigma


!Define E normalizatiopn
Enorm = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2) + delta*(sin(theta)*phi_dot)**2
Eprime2 = (Eobs**2 - (sigma-2.0_dp*r)*omega2/sigma)/Enorm
Eprime = sqrt(Eprime2)




print *, 'IC', r_dot, theta_dot, phi_dot

stop

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





print *, 'transform to global'
print *, '4-vel. Contravar = ', ui
print *, 'IN:', ki
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



print *, u_contra(3), u_covar(3)

delta = r**2.0_dp + a**2.0_dp - 2.0_dp*r
N1 = sqrt(- metric_covar(2,2) * (u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4)) * (1.0_dp + u_covar(3)*u_contra(3)) )
N2 = sqrt(metric_covar(3,3) * (1.0_dp + u_covar(3) * u_contra(3)) )
N3 = sqrt(-(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))*delta*sin(theta)**2)



print *, 'Calc = ', ki(3)/N2

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




print *, 'out:', ki_global

!A check
call mag_4(r,theta,ki_global,mag)
print *, 'Magnitude checks'
print *, -(ki(1))**2 + ki(2)**2 + ki(3)**2 + ki(4)**2
print *, mag
ki = ki_global

end subroutine transform_to_global





subroutine RT_Forward()

        print *, 'for'
end subroutine RT_Forward









subroutine RT_Backward(xi)
!Arguments
real(kind=dp), dimension(4), intent(in) :: xi !position of MPD
!Other
real(kind=dp) :: xT, yT, zT, mm
real(kind=dp) :: alpha,beta
real(kind=dp),dimension(6) :: IO
real(kind=dp),dimension(4) :: globals
real(kind=dp),dimension(3) :: OT
integer(kind=dp) :: plot !Do you want to save ray path?


!The target points 
mm = sqrt(xi(2)**2 + a**2)
xT = mm * sin(xi(3))*cos(xi(4))
yT = mm * sin(xi(3))*sin(xi(4))
zT = mm * cos(xi(3))


xT = -13.0_dp
yT = 5.0_dp
zT = 0.0_dp




!TO DO()
!Secondary rays
!Speed optimization
!Robustness for edge cases e.g. hits BH

!Initial guess at the coordinates of the image plane
alpha = -10.0_dp*yT !Falls below horizon for alpha=-yT (=-5.0_dp)
beta = xT*cos(ThetaObs) + zT*sin(ThetaObs)

!Get the target point and initial alpha beta
IO(1) = xT
IO(2) = yT
IO(3) = zT
IO(4) = alpha
IO(5) = beta 

!Fire a ray. How it does i.e. ds is written to IO(6)
plot = 0
call shoot(IO,plot)

print *, 'Initial ds =', IO(6)

!Now run some conjugate gradient descent (non-linear) using the AIM subroutine
OT(1) = 0.010_dp ! Initial stepsize alpha for backtracing
globals = 0.0_dp

do while (IO(6) .GT. 1.0d-6)

call aim(IO,globals,OT)

enddo


!Trace the ray - useful for data+plotting. Or just write to file for use later

!The mininum has been sucessfully found. Now tracing the ray
print *, 'The mininum has been sucessfully found. Now tracing the ray'
plot=1
call shoot(IO,plot)

stop



end subroutine RT_Backward


subroutine aim(IO,globals,OT)
!Arguments
real(kind=dp),dimension(6),intent(inout) :: IO
real(kind=dp),dimension(4),intent(inout) :: globals
real(kind=dp),dimension(:),intent(inout) :: OT !i.e. other
!Other
real(kind=dp),parameter :: gbit = 1.0e-18_dp
real(kind=dp) :: dsA, dsB, gA, gB, dsORIG, zeta,hA, hB
real(kind=dp) :: etaRT, factor, alpha0, beta0,dsBEST,aBEST,bBEST
real(kind=dp) :: tau, c,t, pA, pB,norm,m, dstemp !Armijo
integer(kind=dp), parameter :: plot = 0 !never want to plot when  aiming.





!Load original ds
alpha0 = IO(4)
beta0 = IO(5)
dsORIG = IO(6)




!Alpha gradient gA
IO(4) = IO(4) + gbit
call shoot(IO,plot)
dsA = IO(6)
gA = (dsA - dsORIG)/gbit
gA = - gA 
IO(4) = IO(4) - gbit



!Beta gradient gB
IO(5) = IO(5) + gbit
call shoot(IO,plot)
dsB = IO(6)
gB = (dsB - dsORIG)/gbit
gB = -gB 
IO(5) = IO(5) - gbit







!Now get the alpha/beta adjustment direction

if (globals(1) .EQ. 0.0_dp) then
!It is the first time, just set zeta = 0
zeta = 0.0_dp
else
zeta = (gA*gA +gB*gB)/(globals(1)**2 + globals(2)**2)
endif


hA = gA +zeta*globals(3) 
hB = gB +zeta*globals(4)


!Got the direction. Now perform a line search

!set inital optimization params
etaRT = OT(1)
factor = 2.0_dp

!loop line search

dsBEST = 1d20
aBEST = alpha0
bBEST = beta0






tau = 0.10_dp
c = 0.0010_dp

!Testing backtrack line search

norm = sqrt(hA**2 + hB**2)
pA = hA/norm
pB = hB/norm

m = gA*pA + gB*pB !Is this the correct defenition?
t = -c*m

etaRT = OT(1)
dstemp = IO(6)


!print *, alpha0, etaRT*hA
!stop

02 do

!if success on first go, try a larger step size?

        IO(4) = alpha0 + etaRT*hA
        IO(5) = beta0 + etaRT*hB
 

        call shoot(IO,plot)

        print *, 'i:', IO(4:6)

       ! if ((IO(6) - dstemp) .LT. etaRT*t) then
        
        if (IO(6) .LT. dstemp) then


        !Found something good. EXIT



!Update globals
globals(1) = gA
globals(2) = gB
globals(3) = hA
globals(4) = hB



if (etaRT .LT. 1e-10_dp) then
!etaRT is very small.
!Reset to normal
OT(1) =1.0_Dp
else
OT(1) = etaRT  !etaRT alwas decas Is this good?
endif




print *, '------'
        return



        else


        etaRT = etaRT*tau

        if (etaRT .LT. 1e-10_dp) then
        !Reset
        IO(6) = dsORIG
        IO(5) = beta0
        IO(4) = alpha0
        globals = 0.0_dp
        OT(1) = 1.0_dp
        return
        endif

        endif





enddo






stop




!01 do
!
!        etaRT =etaRT * factor
!        IO(4) = aBEST + etaRT*hA
!        IO(5) = bBEST + etaRT*hB
!
!
!        call shoot(IO)
!
!
!
!
!        if (IO(6) .LT. dsBEST) then
!        aBEST = IO(4)
!        bBEST = IO(5)
!        dsBEST = IO(6)
!
!        else
!                EXIT !Exit do loop
!        endif
!
!
!
!
!
!
!enddo
!
!
!
!
!
!if (dsBEST .LT. dsORIG) then
!!Good. Update
!
!
!!Update search attempts        
!IO(4) = aBEST
!IO(5) = bBEST
!IO(6) = dsBEST
!
!
!!Update globals
!globals(1) = gA
!globals(2) = gB
!globals(3) = hA
!globals(4) = hB
!
!
!!Update stepsize
!OT(1) = etaRT/5.0_dp
!
!
!
!else
!
!!Dont update
!IO(4) = alpha0
!IO(5) = beta0
!IO(6) = dsORIG
!
!!Set search parameters
!OT(1) = OT(2) !etaRT = etaLOW
!OT(3) = OT(3) + 1 !Track the number of fails
!
!endif
!
!
!
!if (OT(3) .GT. 20) then
!!Repeatedly fails to lower
!!Change search parameters and reset globals
!OT(1) = 5.0e-9_dp
!OT(2) = OT(1)
!OT(3) = 0.0_dp
!globals = 0.0_dp
!endif
!
!

!print *, 'OUT:', IO(6), OT(1), OT(3)


end subroutine aim




subroutine shoot(IO,plot)
!Arguments
real(kind=dp),dimension(6),intent(inout) :: IO
integer(kind=dp), intent(in) :: plot
!Other
real(kind=dp) :: xT, yT, zT, xP, yP, zP,mm
real(kind=dp) :: alpha,beta
real(kind=dp) :: xprime, yprime, zprime
real(kind=dp) :: w, r, theta,phi,t
real(kind=dp) :: rdot, thetadot, phidot, zdot, sig,u,vv
real(kind=dp), dimension(7) :: ray !Ray initial conditions
real(kind=dp), dimension(6) :: v !The variables e.g. r, theta phi etc
real(kind=dp), dimension(4) :: c !The constants L. kappa etc
real(kind=dp), dimension(2) :: b !the BackArray - contains stuff which relates to intersections and backwards ray tracing
integer(kind=dp) :: counter,i
real(kind=dp) :: Rstart, xOUT,dx, ds2
integer(kind=dp), parameter :: nrows = 1e8,ncols = 3
real(kind=dp), dimension(:,:), allocatable :: PlotArray
character(len=300) :: RTFile


if (plot .EQ. 1) then
!Set up to plot
allocate(PlotArray(nrows,ncols))
call create_RT_plotfile(IO(4),IO(5),RTFile)
endif

!Read in array
xT = IO(1)
yT = IO(2)
zT = IO(3)
alpha = IO(4)
beta = IO(5)




!Convert to the primed Cartesian frame
xprime = sqrt(Robs**2.0_dp +a**2.0_dp) * sin(ThetaObs) - beta*cos(ThetaObs)
yprime = alpha
zprime = Robs*cos(ThetaObs) + beta*sin(ThetaObs)




!Convert it to Boyer Lindquist
w = xprime**2.0_dp +yprime**2.0_dp +zprime **2.0_dp - a**2.0_dp
r = sqrt((w+sqrt(w**2.0_dp+4.0_dp*a**2.0_dp*zprime**2.0_dp))/2.0_dp)
theta = acos(zprime/r)
phi = atan2(yprime,xprime)
t=0.0_dp


!And then get the derivatives
sig = r**2.0_dp +(a*cos(theta))**2.0_dp
u = sqrt(r**2.0_dp+a**2.0_dp)
vv= -sin(ThetaObs)*cos(phi)
zdot = -1.0_dp



rdot = -zdot*(-u**2.0*cos(ThetaObs)*cos(theta)+r*u*vv*sin(theta))/sig
thetadot = -zdot*(cos(ThetaObs)*r*sin(theta)+u*vv*cos(theta))/sig	
phidot = -zdot*sin(ThetaObs)*sin(phi)/(u*sin(theta))



!Write to array
ray(1) = t
ray(2) = r
ray(3) = theta
ray(4) = phi
ray(5) = rdot
ray(6) = thetadot
ray(7) = phidot

!Set up initial conditions
call GeneralInitialConditions(ray,v,c)
Rstart = v(1)
counter = 1

b(1) = xT
b(2) = 0.0_dp







!Iterate
do while (v(1) .GT. Rhor .and. v(1) .LT. rstart*10.0_dp)
!Do a timestep

call rk_geodesic(v,c,b)
!print *, plot, v(1)
if (plot .EQ. 1) then
!Save to array for plotting
PlotArray(counter,1:3) = v(1:3) 
endif
counter = counter + 1


!Check to see if condition is satisfied
if (b(2) .EQ. 1.0_dp) then

!Calculate ds2
mm = sqrt(v(1)**2 + a**2)
xP = mm*sin(v(2))*cos(v(3))
yP = mm*sin(v(2))*sin(v(3))
zP = v(1)*cos(v(2))

ds2 = (xP - xT)**2 + (yP - yT)**2 + (zP-zT)**2
IO(6) = ds2



!Save the plot array


if (plot .EQ. 1) then

!MAKE THIS A SUBROUTINE

open(unit =10, file = RTFile,form='formatted')
do i=1,counter-1
mm = sqrt(PlotArray(i,1)**2 + a**2)
xP = mm*sin(PlotArray(i,2))*cos(PlotArray(i,3))
yP = mm*sin(PlotArray(i,2))*sin(PlotArray(i,3))
zP = PlotArray(i,1)*cos(PlotArray(i,2))
write(10,*) xP , yP, zP
enddo
close(10)
print *, 'saved'
deallocate(PlotArray)
endif





return
endif

enddo



print *, 'Fell below horizon'
stop






!open(unit =10, file = RTFile,form='formatted')
!do i=1,counter-1
!mm = sqrt(PlotArray(i,1)**2 + a**2)
!xP = mm*sin(PlotArray(i,2))*cos(PlotArray(i,3))
!yP = mm*sin(PlotArray(i,2))*sin(PlotArray(i,3))
!zP = PlotArray(i,1)*cos(PlotArray(i,2))
!write(10,*) xP , yP, zP
!enddo
!close(10)
!print *, 'saved'
!deallocate(PlotArray)

end subroutine shoot





subroutine GeneralInitialConditions(ray,v,c)
!Argument
real(kind=dp), dimension(7) :: ray !Ray initial conditions
real(kind=dp), dimension(6), intent(out) :: v
real(kind=dp), dimension(4), intent(out) :: c
!Other
real(kind=dp) :: r,theta,rdot,thetadot,phidot, sigma,delta
real(kind=dp) :: E2, E, Enorm, Eprime2, Eprime,Eobs
real(kind=dp) :: B2, fr,ft,omega2
real(kind=dp) :: Lz, pr, ptheta, phi,kappa, En, s1

!Load the data
r = ray(2)
theta = ray(3)
phi = ray(4)
rdot = ray(5)
thetadot = ray(6)
phidot = ray(7)

!Setup some defenitions
sigma = r**2 +a**2*cos(theta)**2
delta = r**2 -2.0_dp*r + a**2



!Define the plasma frequency
B2 =  N*4.0_dp*PI*electron_charge**2 / electron_mass
call plasma_fr(r,fr)
call plasma_ft(r,ft)
omega2 = B2 * (fr+ft)/sigma

!Construct the energy is the conserved energy?
!Eobs = 16.0_dp
!Enorm = (sigma-2.0_dp*r)*(rdot**2/delta + thetadot**2) + delta*(sin(theta)*phidot)**2
!Eprime2 = (Eobs**2)/Enorm
!Eprime = sqrt(Eprime2)

!Correct to ensure correct Energies
!rdot = rdot * Eprime
!thetadot = thetadot*Eprime
!phidot = phidot*Eprime


!If you want you can check that E = Eobs
!E2 = (sigma-2.0_dp*r)*(rdot**2/delta + thetadot**2) + delta*(sin(theta)*phidot)**2
!E = sqrt(E2)


! Compute the energy
s1 = sigma-2.0*r
E2 = s1*(rdot**2.0/delta +thetadot**2.0 +omega2/sigma) + delta*sin(theta)**2.0*phidot**2.0
En  = sqrt(E2)






!Get the angular momentum
Lz = (sigma*delta*phidot - 2.0_dp*a*r*En)*sin(theta)**2 / (s1)


!Get the momenta (non constant)
pr = rdot * sigma/delta
ptheta = sigma*thetadot

!Normalise to E = 1
pr = pr/En
ptheta = ptheta/En
Lz = Lz/En
B2 = B2 / E2




!Define one last constant and export
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
c(4) = 1.0d-3 !Initial stepsize for RT


!print *, 'Iliad IC:'
!print *, v(1:3)
!print *, v(4:6)
!print *, c, En



end subroutine GeneralInitialConditions






end module RayTracing
