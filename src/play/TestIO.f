program TEST_IO


implicit none

!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)

!Define array dimensions
integer(kind=dp) :: nrows = 1d4, ncols = 12

!Define write/read arrays
real(kind=dp), dimension(:,:), allocatable :: ArrayWrite, ArrayRead, ArrayRead2


!Allocate
allocate(ArrayWrite(nrows,ncols))
allocate(ArrayRead(nrows,ncols))
allocate(ArrayRead2(nrows,ncols))

!Populate the array
ArrayWrite=1.0_dp

!Write to file

open(unit=10, file='example.dat',status='replace', form='unformatted',access='stream')
write(10) ArrayWrite
close(10)


!Re-populate array
ArrayWrite=2.0_dp

!And append to existing file
open(unit=10, file='example.dat', form='unformatted',position='append',access='stream')
write(10) ArrayWrite
close(10)

!Now read in all the data

open(unit=10, file='example.dat' , form='unformatted',access='stream')
read(10) ArrayRead
read(10) ArrayRead2
close(10)
print *, ArrayRead(nrows,1)
print *, ArrayRead2(nrows,1)


end program TEST_IO
