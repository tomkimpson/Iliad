program main

use parameters
use OrbitalDynamics
use RayTracing
implicit none



!Set up some paths, save directories, etc

call setup()



!Generate the orbit
!print *, 'Calling MPD'
!call MPD()




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

character(len = 20), parameter :: FMT1 = "(F10.2)"
character(len = 20), parameter :: FMT3 = "(F10.3)"
character(len = 80) :: BHMass, BHSpin, PSRMass, PSRSpin, PSRangle1, PSRangle2, PSRangle3, OrbitalPeriod, OrbitalEcc, &
                       OrbitalNumber, OrbitalIota
call get_environment_variable("IliadDir", path)
MPDBinaryData = trim(adjustl(path))//'MPDBinaryData.dat'
MPDFormatData = trim(adjustl(path))//'MPDFormatData.txt'
RTPath = trim(adjustl(path))//'RT/'

!Welcome messages

print *, 'Iliad is running'

!Write to string
write(BHMass, FMT1) MBH/1.0e6_dp
write(BHSpin, FMT1) a

write(PSRMass, FMT1) MPSR
write(PSRSpin, FMT1) p0*1.0e3_dp

write(PSRangle1, FMT1) stheta
write(PSRangle2, FMT1) sphi
write(PSRangle3, FMT1) chi

write(OrbitalPeriod, FMT3) KeplerianPeriod
write(OrbitalEcc, FMT1) eccentricity
write(OrbitalNumber, FMT1) Norbit
write(OrbitalIota, FMT1) iota


print *, 'You have selected the following settings:'

print *, '-------------------------'
print *, '-------------------------'
print *, 'BH PARAMETERS:'
print *, '-------------------------'
print *, '-------------------------'
print *, 'BH mass =', trim(adjustl(BHMass)), ' x 10^6 solar masses. BH spin = ', trim(adjustl(BHSpin))
print *, '-------------------------'
print *, '-------------------------'
print *, 'PSR PARAMETERS:'
print *, '-------------------------'
print *, '-------------------------'
print *, 'PSR mass =', trim(adjustl(PSRMass)), ' solar masses. PSR spin period = ', trim(adjustl(PSRSpin)), ' milliseconds'
print *, 'Stheta = ', trim(adjustl(PSRangle1)), ' Sphi = ', trim(adjustl(PSRangle2)), 'Chi =', trim(adjustl(PSRangle3)) 
print *, '-------------------------'
print *, '-------------------------'
print *, 'ORBITAL PARAMETERS:'
print *, '-------------------------'
print *, '-------------------------'
print *, 'Period = ', trim(adjustl(OrbitalPeriod)), ' years. Eccentricity = ', trim(adjustl(OrbitalEcc))
print *, 'Inclination = ', trim(adjustl(OrbitalIota)), ' Number of orbits = ', trim(adjustl(OrbitalNumber))


if (lambda .EQ. 1) then
print *, 'Spin-curvature coupling is turned on'
else
print *, 'Spin-curvature coupling is turned off'
endif

print *, '---- END INITIAL SETUP ----'
end subroutine setup


