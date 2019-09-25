module IO

use parameters
use constants



implicit none

private

public FileOpen, ToTextFile

contains


subroutine ToTextFile(f)

character(len=200) :: f
real(kind=dp), dimension(13) :: row
real(kind=dp), dimension(:,:), allocatable :: array
integer(kind=dp) :: i,stat


print *, 'Called ToTExtFile'
allocate(array(nrows,entries+1))

!Open the .dat file
open(unit=10, file=MPDBinaryData , form='unformatted',access='stream')

!Open the .txt file
open(unit=20, file=MPDFormatData , form='formatted',access='stream')




stat = 0

do while (stat .EQ. 0)

    read(10,iostat=stat) array
    
    print *, array, nrows
    do i=1,nrows
    
    if (array(i,2) .EQ. 0.0_dp) then
    exit
    endif

    print *, 'Writing', MPDFormatData
    
    write(20,*) array(i,2), array(i,3), array(i,4), a
    
    
    enddo
    array = 0.0_dp



enddo


close(10)
close(20)




stop



end subroutine ToTextFile


subroutine FileOpen(f)
character(len=200) :: f
logical :: res


inquire( file=f, exist=res )
if (res) then

open(unit=10, file=MPDBinaryData,status='old',position='append', form='unformatted',access='stream')
        

else

open(unit=10, file=MPDBinaryData,status='replace', form='unformatted',access='stream')


endif


end subroutine FileOpen




end module IO
