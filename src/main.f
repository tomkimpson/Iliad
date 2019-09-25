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
stop


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




end subroutine setup


