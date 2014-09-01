program analyse
use TrajReader
implicit none

!needed to open the trajectory file
logical :: file_opened =  .False.
character(len=128) :: filename
integer :: filename_len

!data to be read from the trajectory file
integer(c_int) timestep, N
real(c_double), dimension(3) :: box_lo, box_hi, prd
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
real(c_double), parameter :: cutoff = 3.0
real(c_double), parameter :: cutoffsq = cutoff*cutoff

integer :: timesteps_read = 0

!interparticle distance
real(c_double), dimension(3) :: del
real(c_double) :: rsq

!intermediate variables
real(c_double) :: tconF_numerator = 0.0
real(c_double) :: tconF_denominator = 0.0
real(c_double) :: tcon1 = 0.0
real(c_double) :: tkin = 0.0
real(c_double) :: tconF = 0.0

!accumulators for averages over entire simulation
real(c_double) :: tkin_ave = 0.0
real(c_double) :: tconF_numerator_ave = 0.0
real(c_double) :: tconF_denominator_ave = 0.0
real(c_double) :: tcon1_ave = 0.0
real(c_double) :: tkinsq_ave = 0.0
real(c_double) :: tconF_numerator_sq_ave = 0.0
real(c_double) :: tconF_denominator_sq_ave = 0.0
real(c_double) :: tcon1sq_ave = 0.0

!error variables
real(c_double) :: tkin_err = 0.0
real(c_double) :: tconF_numerator_err = 0.0
real(c_double) :: tconF_denominator_err = 0.0
real(c_double) :: tconF_err = 0.0
real(c_double) :: tcon1_err = 0.0

call get_command_argument(1, filename, filename_len)

file_opened = OpenFile(filename, filename_len)

if (.not. file_opened) then
        print*, "Couldn't open file ", filename
else
        open(unit=output_file, file=output_filename, action="write", status="replace")
        write(output_file,*) "# Timestep   Tkin   Tcon"
        do
                loop_continue = readframe(timestep, N, box_lo, box_hi, x, v, f)
                if (.not. loop_continue) exit
                print*, timestep

                ! calculate the periodic replica distance
                ! needed for the minimum_image subroutune
                prd = box_hi - box_lo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!! START MODIFYING HERE !!!

                ! the following have been read from the trajectory for you:
                ! integer timestep
                ! integer N (number of atoms)
                ! real box_lo(3) minimum coordinates in all three directions
                ! real box_hi(3) maximum coordinates in all three directions
                ! real prd(3) periodic replica distance in all three directions
                ! real x(N,3) atomic positons
                ! real v(N,3) atomic velocities
                ! real f(N,3) total force on each atom F_i

                ! the following subprogrammes are provided:
                ! - minimum_image(vector, prd) : apply the minimum image convention
                ! - veclengthsq(vector) : return the squared length of the vector x**2 + y**2 + z**2

                ! You need to give values to the following variables:
                ! - tkin
                ! - tconF_denominator
                ! - tconF_numerator
                ! the existing code will take care of the averages and standard deviations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do i=1,N
                        tkin = tkin + veclengthsq(v(i,:))
                        tconF_denominator = tconF_denominator + veclengthsq(f(i,:))
                        do j=i+1,N
                                del = x(i,:) - x(j,:)
                                call minimum_image(del, prd) 
                                rsq = veclengthsq(del)

                                if (rsq < cutoffsq) then
                                        ! we add the contributions for i->j and j->i
                                        tconF_numerator = tconF_numerator - 2.0*div_f(rsq)
                                end if
                        end do
                end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!! STOP MODIFYING HERE !!!
                ! (Unless you wish to add a function/subroutine) !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                tkin = mass * tkin / (3.0 * N)
                tkin_ave = tkin_ave + tkin
                tkinsq_ave = tkinsq_ave + tkin*tkin

                tcon1 = tconF_denominator / tconF_numerator
                tcon1_ave = tcon1_ave + tcon1
                tcon1sq_ave = tcon1sq_ave + tcon1*tcon1

                tconF_numerator_ave = tconF_numerator_ave + tconF_numerator
                tconF_numerator_sq_ave = tconF_numerator_sq_ave + tconF_numerator*tconF_numerator
                tconF_denominator_ave = tconF_denominator_ave + tconF_denominator
                tconF_denominator_sq_ave = tconF_denominator_sq_ave + tconF_denominator*tconF_denominator
                !write the current temperatures to the output file
                write(output_file,*) timestep, tkin, tcon1

                timesteps_read = timesteps_read + 1

                !reset the instantaneous variables for next time
                tconF_numerator = 0.0
                tconF_denominator = 0.0
                tcon1 = 0.0
                tkin = 0.0
                tconF = 0.0
        end do
        close(output_file)

        tkin = tkin_ave / timesteps_read
        tkin_err = stderr(tkin, tkinsq_ave / timesteps_read, timesteps_read)
        tcon1 = tcon1_ave / timesteps_read
        tcon1_err = stderr(tcon1, tcon1sq_ave / timesteps_read, timesteps_read)
        tconF_numerator = tconF_numerator_ave / timesteps_read
        tconF_numerator_err = stderr(tconF_numerator, tconF_numerator_sq_ave / timesteps_read, timesteps_read)
        tconF_denominator = tconF_denominator_ave / timesteps_read
        tconF_denominator_err = stderr(tconF_denominator, tconF_denominator_sq_ave / timesteps_read, timesteps_read)

        tconF = tconF_denominator / tconF_numerator
        tconF_err = tconF**2 * ((tconF_numerator_err/tconF_numerator)**2 + (tconF_denominator_err/tconF_denominator)**2)

        print*, "Read", timesteps_read, "timesteps."
        print*, "Temperatures"
        print*, "Kinetic: ", tkin, "+/-", tkin_err 
        print*, "Con1", tcon1, "+/-", tcon1_err
        print*, "ConF", tconF, "+/-", tconF_err
end if
contains
real(c_double) function div_f(rsq)
        real(c_double), intent(in) :: rsq
        real(c_double) :: r6, r8, r14

        r6 = rsq*rsq*rsq
        r8 = r6*rsq
        r14 = r8*r6

        div_f = 30.0*C6/r8 - 132.0*C12/r14
end function div_f
real(c_double) function stderr(X, X2, N)
        implicit none
        real(c_double), intent(in) :: X, X2
        integer(c_int), intent(in) :: N

        stderr = X2 - X**2
        stderr = sqrt(stderr/N)
end function stderr
real(c_double) function veclengthsq(vec)
        implicit none
        real(c_double), dimension(3), intent(in) :: vec

        veclengthsq = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
end function veclengthsq
subroutine minimum_image(vec, prd)
        real(c_double), dimension(3), intent(inout) :: vec
        real(c_double), dimension(3), intent(in) :: prd
        integer i

        do i=1,3
                if (vec(i) < -prd(i)/2.0) then
                        vec(i) = vec(i) + prd(i)
                else if (vec(i) > prd(i)/2.0) then
                        vec(i) = vec(i) - prd(i)
                end if
        end do
end subroutine minimum_image
end program analyse
