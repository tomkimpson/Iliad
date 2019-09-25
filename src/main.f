program main

use parameters
use OrbitalDynamics
use RayTracing
implicit none



!Set up some paths, save directories, etc

call setup()



!Generate the orbit
print *, 'Calling MPD'
call MPD()


!Use the orbit as initial conditions for ray tracing. 
print *, 'Calling RayTracing'
call RT()





print *, 'Completed'
print *, 'Output files are:'
print *, MPDFormatData
print *, RTPath


end program main






subroutine setup()

use constants

call get_environment_variable("IliadDir", path)
MPDBinaryData = trim(adjustl(path))//'MPDBinaryData.dat'
MPDFormatData = trim(adjustl(path))//'MPDFormatData.txt'
RTPath = trim(adjustl(path))//'RT/'


!Welcome messages

print *, 'Iliad is running'

print *, 'You have selected the following settings:'

print *, '-------------------------'
print *, 'BH:'
print *, '-------------------------'
print *, 'BH mass =', MBH/1.0d6, ' BH spin = ', a
print *, '-------------------------'
print *, 'PSR:'
print *, '-------------------------'
print *, 'PSR mass =', MPSR, ' PSR spin period = ', p0
print *, 'Stheta = ', stheta, 'Sphi = ', sphi
print *, 'Chi angle = ', chi
print *, '-------------------------'
print *, 'Orbit:'
print *, '-------------------------'
print *, 'Period = ', KeplerianPeriod, ' years. Eccentricity = ', e
print *, 'Inclination = ', iota, ' Number of orbits = ', Norbit

if (lambda .EQ. 1) then
print *, 'Spin-curvature coupling is turned on'
else
print *, 'Spin-curvature coupling is turned off'
endif



end subroutine setup


