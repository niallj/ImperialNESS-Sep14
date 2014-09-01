      SUBROUTINE READ_INPUTDATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Reads in input data from a file
!!!!  uses a namelist
!!!!  Karl Travis 17th Jan 2012 (fortran 90/95 version) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      USE header_file
      IMPLICIT NONE

      OPEN(3,FILE = 'input.dat')

      READ(UNIT = 3, NML = input_deck)

      WRITE(UNIT = *, NML = input_deck) ! echo the stuff to screen for user to check

      CLOSE(3)

      SELECT CASE (IPOT_TYPE) ! choose option depending on potential type

        CASE(1)   ! Cut LJ potential

          CUTSQ = RCUT**2

        CASE(2)   ! Splined-LJ potential

          RCUT = 1.737051855D00 
          CUTSQ = RCUT**2
          RCUT_INNER = 1.24445506D0 
          RCUT_INNER_SQ = RCUT_INNER**2 
          SPLINE_R2  = -4.864890083D00 
          SPLINE_R3  = -3.2920028D00 

        CASE(3)   ! Smooth repulsive potential

          RCUT = SIGMA
          CUTSQ = RCUT**2
          EPSIG8 = EPS*8.0D0/CUTSQ

        CASE(4)   ! m-n family of potentials

          IF(M_POWER.GE.N_POWER) THEN
            WRITE(*,*) 'Oops: cant have m >= n. Try again!'
            STOP 
          END IF
          DIFF_MN = N_POWER - M_POWER
          M_RATIO = DBLE(M_POWER)/DIFF_MN
          N_RATIO = DBLE(N_POWER)/DIFF_MN
          FORCE_MN_FACT = 2.0D0*DBLE(M_POWER*N_POWER)/DIFF_MN
          N_POWER_SUB1 = N_POWER - 1
          M_POWER_SUB1 = M_POWER - 1
          CUTSQ = 2.0d0
          RCUT = SQRT(CUTSQ)

        CASE(5)   ! Lucy potential

          CUTSQ = RCUT**2
          EFACT_LUCY = 5.0D0/(CUTSQ*PI)
          FORCE_FACT_LUCY = EFACT_LUCY*12.0D0/CUTSQ

        CASE(6)   ! Soft disk potential

          SOFT_N = DBLE(SS_POWER)
          CUTSQ = RCUT**2

        CASE(7)   ! WCA potential 
          RCUT = 2**(1.0d0/6.0d0)
          CUTSQ = RCUT**2

        CASE DEFAULT 

          WRITE(*,*) 'Supports only IPOT_TYPE values from 1..7 at the moment'
          STOP

      END SELECT
       TAU_TEMP_SQD = TAU_TEMP**2

      END SUBROUTINE
