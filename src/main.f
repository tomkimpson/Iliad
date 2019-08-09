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
end program main






subroutine setup()

use constants

call get_environment_variable("IliadDir", path)

MPDBinaryData = trim(adjustl(path))//'MPDBinaryData.dat'
MPDFormatData = trim(adjustl(path))//'MPDFormatData.txt'

end subroutine setup


