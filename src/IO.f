module IO

use parameters
use constants



implicit none

private

public FileOpen

contains

subroutine ReadFile(f)

character(len=200) :: f
integer(kind=dp) :: i,stat,stat_io, j
real(kind=dp), dimension(:,:), allocatable :: ReadArray

allocate(ReadArray(nrows,entries))

open(unit = 10,FORM = 'unformatted', file = f, iostat= stat,STATUS = 'OLD', access='stream')

read(10) ReadArray
read(10) ReadArray
close(10)

end subroutine ReadFile


subroutine FileOpen(f)
character(len=200) :: f
logical :: res


inquire( file=f, exist=res )
if (res) then
print *, 'File Exists'

open(10,file=f,status="old",action="write",form='unformatted',position='append',access='stream')


else
print *, 'new file'
open(unit = 10,file=f,status="new",action="write",form='unformatted',access='stream')
endif


end subroutine FileOpen




end module IO
