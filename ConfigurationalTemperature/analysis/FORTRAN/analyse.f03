program analyse
use TrajReader
implicit none

!needed to open the trajectory file
logical :: file_opened =  .False.
character(len=128) :: filename
integer :: filename_len

!data to be read from the trajectory file
integer(c_int) timestep, N
real(c_double), dimension(3) :: box_lo, box_hi
real(c_double), dimension(:,:), allocatable :: x, v, f
logical(c_bool) :: loop_continue = .True.

!file output control
integer, parameter :: output_file = 1
character(*), parameter :: output_filename = "temperatures.dat"

!loop counters
integer :: i, j

!values for temperature calculations
real(c_double), parameter :: mass = 1.0
real(c_double), parameter :: C6 = 4.0
real(c_double), parameter :: C12 = 4.0

integer :: timesteps_read = 0

!intermediate variables
real(c_double) :: tcon_numerator = 0.0
real(c_double) :: tcon_denominator = 0.0
real(c_double) :: tconF = 0.0
real(c_double) :: tcon1 = 0.0
real(c_double) :: tkin = 0.0

!accumulators for averages over entire simulation
real(c_double) :: tkin_ave = 0.0
real(c_double) :: tconF_numerator_ave = 0.0
real(c_double) :: tconF_denominator_ave = 0.0
real(c_double) :: tcon1_ave = 0.0
real(c_double) :: tkinsq_ave = 0.0
real(c_double) :: tconF_numerator_sq_ave = 0.0
real(c_double) :: tconF_denominator_sq_ave = 0.0
real(c_double) :: tcon1sq_ave = 0.0

call get_command_argument(1, filename, filename_len)

file_opened = OpenFile(filename, filename_len)

if (.not. file_opened) then
        print*, "Couldn't open file ", filename
else
        open(unit=output_file, file=output_filename, action="write", status="replace")
        write(output_file,*) "# Timestep   Tkin   Tcon1   TconF"
        do
                loop_continue = readframe(timestep, N, box_lo, box_hi, x, v, f)
                if (.not. loop_continue) exit

                !the following have been read from the trajectory for you:
                ! integer timestep
                ! integer N (number of atoms)
                ! real box_lo(3) minimum coordinates in all three directions
                ! real box_hi(3) maximum coordinates in all three directions
                ! real x(N,3) atomic positons
                ! real v(N,3) atomic velocities
                ! real f(N,3) total force on each atom F_i

                !remember that you need to calculate:
                ! tkin = the standard "equipartition" temperature
                ! tcon1, tconF, as defined in the booklet
                ! we need instantaneous values (instantaneously, tcon1 = tconF)
                ! and also averages over the whole simulation
                ! for the final averages, you should also calculate an
                ! error, using Var(X) = <X**2> - <X>**2, and
                ! stderr = sqrt(Var(X)/N)
 

                !write the current temperatures to the output file
                write(output_file,*) timestep, tkin, tconF, tcon1

                timesteps_read = timesteps_read + 1
        end do
        close(output_file)
        print*, "Read", timesteps_read, "timesteps."
end if
end program analyse
