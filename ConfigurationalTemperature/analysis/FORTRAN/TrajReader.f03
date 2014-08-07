module TrajReader
use, intrinsic :: iso_c_binding
implicit none
interface
logical(c_bool) function FrameBody(n, x, v, f) bind(c)
        use, intrinsic :: iso_c_binding
        implicit none
        real(c_double), dimension(:,:), intent(out) :: x, v, f
        integer(c_int) :: n
end function FrameBody

logical(c_bool) function FrameHeader(timestep, n, boxlo, boxhi) bind(c)
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), intent(out) :: timestep, n
        real(c_double), dimension(3), intent(out) :: boxlo, boxhi
end function FrameHeader

logical(c_bool) function OpenFile(filename, len) bind(c)
        use, intrinsic :: iso_c_binding
        implicit none
        character(c_char), dimension(*):: filename
        integer(c_int) :: len
end function OpenFile
end interface

contains
logical(c_bool) function readframe(timestep, n, boxlo, boxhi, x, v, f)
        implicit none
        integer(c_int), intent(out) :: timestep, n
        real(c_double), dimension(3), intent(out) :: boxlo, boxhi
        real(c_double), dimension(:,:), allocatable, intent(out) :: x, v, f

        !first clean up the arrays
        if (allocated(x)) then
                deallocate(x)
        end if
        if (allocated(v)) then
                deallocate(v)
        end if
        if (allocated(f)) then
                deallocate(f)
        end if

        readframe = FrameHeader(timestep, n, boxlo, boxhi)

        if (readframe) then
                ! allocate the arrays
                allocate(x(n, 3))
                allocate(v(n, 3))
                allocate(f(n, 3))
                readframe = FrameBody(n, x, v, f)
        end if
end function readframe
end module TrajReader
