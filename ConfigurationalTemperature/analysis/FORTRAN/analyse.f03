program analyse
use TrajReader
implicit none

logical :: file_opened =  .False.
character(len=128) :: filename
integer :: filename_len

integer(c_int) timestep, N
real(c_double), dimension(3) :: box_lo, box_hi
real(c_double), dimension(:,:), allocatable :: x, v, f
logical(c_bool) :: loop_continue = .True.

call get_command_argument(1, filename, filename_len)

file_opened = OpenFile(filename, filename_len)

if (.not. file_opened) then
        print*, "Couldn't open file ", filename
else
        do
                loop_continue = readframe(timestep, N, box_lo, box_hi, x, v, f)
                if (.not. loop_continue) exit
                print*, timestep
                end do
end if
end program analyse
