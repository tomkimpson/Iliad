module NumericalMethods


use parameters
use constants
use derivatives


implicit none

private adaptive_shrink, adaptive_grow

public rk, rk_geodesic
contains


subroutine rk_geodesic(y,c,b)
!Arguments
! http://astroa.physics.metu.edu.tr/MANUALS/intel_ifc/mergedProjects/optaps_for/fortran/optaps_prg_arrs_f.htm
real(kind=dp), dimension(:) :: y !Input vector of general size. 
real(kind=dp), dimension(:) :: c !Some constants and stepsizes. 
real(kind=dp), dimension(:) :: b ! Extras for backwards RT 
!Other
real(kind=dp), dimension(size(y)) :: y1,y2,y3,y4,y5,y6
real(kind=dp), dimension(size(y)) :: dy1,dy2,dy3,dy4,dy5,dy6
real(kind=dp), dimension(size(y)) :: k1,k2,k3,k4,k5,k6
real(kind=dp), dimension(size(y)) :: ynew, yerr,deltaErr, yscal, ratio
real(kind=dp) :: h, errmax,xOUT, xT, rT,dx


!Define some stuff from arrays
xT = b(1) ! X coordinate of the target point
rT = b(2) ! R coordinate. !Dont actually use this!


!Set the stepsize
h = c(4)

11 continue


! Y1
y1 = y
call geodesic(y1,c,dy1)
k1 = h * dy1



!Y2
y2 = y1 + B21*k1
call geodesic(y2,c,dy2)
k2 = h * dy2



!Y3
y3 = y1 + B31*k1 + B32*k2
call geodesic(y3,c,dy3)
k3 = h * dy3



!stop

!Y4
y4 = y1 + B41*k1 + B42*k2 + B43*k3
call geodesic(y4,c,dy4)
k4 = h * dy4


!Y5
y5 = y1 + B51*k1 + B52*k2 + B53*k3 + B54*k4 
call geodesic(y5,c,dy5)
k5 = h * dy5


!Y6
y6 = y1 + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
call geodesic(y6,c,dy6)
k6 = h * dy6

!Update
ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6

deltaErr = abs(ynew - yerr)




yscal = abs(y1) + abs(k1) + 1.0d-3
ratio = deltaErr/yscal
errmax = escal * maxval(ratio)






if (RayTracingDirection .EQ. -1.0_dp) then

    !Backwards RT
        xOUT = sqrt(ynew(1)**2 + a**2)*sin(ynew(2))*cos(ynew(3))
        dx = xOUT - xT


        if (abs(dx) .LT. dx_eps) then

 !       print *, 'Good!', dx, xT

        b(2) = 1.0_dp
        y = ynew
        !Leave this subroutine
        return 
        endif



        if (dx .LT. 0.0_dp) then
        h = h*0.50_dp
        goto 11
        endif


endif







if (errmax .GT. 1) then
!The error is too big. Reduce the step size and exit without updating the variable vector
call adaptive_shrink(errmax,h)
goto 11 !?
else
!The error is OK. Grow the stepsize a little and set the variables for the next integration step

tau = tau + h !Update the proper time for MPD
call adaptive_grow(errmax,h)
y = ynew

endif

c(4) = h !Update stepsize. Only if using adaptive stepsize









end subroutine rk_geodesic






subroutine rk(y,c)
!Arguments
! http://astroa.physics.metu.edu.tr/MANUALS/intel_ifc/mergedProjects/optaps_for/fortran/optaps_prg_arrs_f.htm
real(kind=dp), dimension(:) :: y !Input vector of general size. 
real(kind=dp), dimension(:) :: c !Some constants and stepsizes. 
!Other
real(kind=dp), dimension(size(y)) :: y1,y2,y3,y4,y5,y6
real(kind=dp), dimension(size(y)) :: dy1,dy2,dy3,dy4,dy5,dy6
real(kind=dp), dimension(size(y)) :: k1,k2,k3,k4,k5,k6
real(kind=dp), dimension(size(y)) :: ynew, yerr,deltaErr, yscal, ratio
real(kind=dp) :: h, errmax


11 continue

!Set the stepsize
h = c(4)

! Y1
y1 = y
call derivs(y1,dy1)
k1 = h * dy1



!Y2
y2 = y1 + B21*k1
call derivs(y2,dy2)
k2 = h * dy2



!Y3
y3 = y1 + B31*k1 + B32*k2
call derivs(y3,dy3)
k3 = h * dy3


!Y4
y4 = y1 + B41*k1 + B42*k2 + B43*k3
call derivs(y4,dy4)
k4 = h * dy4


!Y5
y5 = y1 + B51*k1 + B52*k2 + B53*k3 + B54*k4 
call derivs(y5,dy5)
k5 = h * dy5


!Y6
y6 = y1 + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5
call derivs(y6,dy6)
k6 = h * dy6

!Update
ynew = y1 + c1*k1  + c3*k3 + c4*k4  +c6*k6 
yerr = y1 + cbar1*k1 + cbar3*k3 + cbar4*k4 + cbar5*k5 + cbar6*k6


deltaErr = abs(ynew - yerr)
yscal = abs(y1) + abs(k1) + 1.0d-3
ratio = deltaErr/yscal
errmax = escal * maxval(ratio)




if (adaptive .EQ. 1) then



if (errmax .GT. 1) then
!The error is too big. Reduce the step size and exit without updating the variable vector



call adaptive_shrink(errmax,h)

c(4) = h !Update stepsize. Only if using adaptive stepsize

goto 11

else
!The error is OK. Grow the stepsize a little and set the variables for the next integration step

tau = tau + h !Update the proper time for MPD
call adaptive_grow(errmax,h)
y = ynew

c(4) = h !Update stepsize. Only if using adaptive stepsize
endif




endif




end subroutine rk




subroutine adaptive_shrink(errmax,dh)
real(kind=dp) :: errmax, htemp,dh


htemp = S*dh*(errmax**PSHRINK)

dh = sign(max(abs(htemp), 0.10_dp*abs(dh)),dh)

end subroutine adaptive_shrink

subroutine adaptive_grow(errmax,dh)
real(kind=dp) errmax,ERRCON,hnext,dh

ERRCON = (5.0_dp/S)**(1.0_dp/PGROW)

!print *, errmax ,ERRCON
if (errmax .GT. ERRCON) then
    hnext = S*dh*errmax**PGROW
 !   print *, 'errcon', S*errmax**PGROW
else
    hnext = 5.0_dp*dh
endif

dh = hnext
END SUBROUTINE adaptive_grow








end module NumericalMethods

