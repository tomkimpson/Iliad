program main

use parameters
use OrbitalDynamics

implicit none



!Set up some paths, save directories, etc

call setup()



!Generate the orbit
print *, 'Calling MPD'
call MPD()



!Use the orbit as initial conditions for ray tracing. 



end program main






subroutine setup()

use constants

call get_environment_variable("IliadDir", path)

MPDBinaryData = trim(adjustl(path))//'MPDBinaryData.dat'

end subroutine setup


