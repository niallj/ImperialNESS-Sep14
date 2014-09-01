      SUBROUTINE PAIR_FORCE(RSQ,FORCEIJ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculates the force between a pair of atoms 
!! Currently supports 6 different forcelaws 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: FORCEIJ,RSQ 
      REAL(kind = double) :: RRR,RR6,RR12 
      REAL(kind = double) :: R,RDIF,RDIF2,RDIF3 
      REAL(kind = double) :: RATIO
      REAL(kind = double) :: FACTOR

      SELECT CASE (IPOT_TYPE) ! choose option depending on potential type

        CASE(1)   ! Cut LJ potential

          RRR = 1.0D0/RSQ
          RR6 = RRR*RRR*RRR
          RR12 = RR6*RR6
          FORCEIJ = 48.0d0*EPS*(RR12 - 0.5d0*RR6)*RRR 

        CASE(2)   ! Splined-LJ potential

          IF(RSQ.LE.RCUT_INNER_SQ) THEN
            RRR = 1.0D0/RSQ
            RR6 = RRR*RRR*RRR
            RR12 = RR6*RR6
            FORCEIJ = 48.0d0*EPS*(RR12 - 0.5d0*RR6)*RRR 
          ELSE
            R = SQRT(RSQ)
            RDIF = R - RCUT
            RDIF2 = RDIF*RDIF
            RDIF3 = RDIF2*RDIF
            FORCEIJ = -(2.0d0*SPLINE_R2*RDIF + 3.0d0*SPLINE_R3*RDIF2)/R
          END IF

        CASE(3)   ! Smooth repulsive potential

          RATIO = RSQ/CUTSQ
          FORCEIJ = EPSIG8*(1.0D0 - RATIO)**3

        CASE(4)   ! m-n family of potentials

          FACTOR = 2.0D0 - RSQ
          FORCEIJ = FORCE_MN_FACT*(FACTOR**N_POWER_SUB1 - FACTOR**M_POWER_SUB1)

        CASE(5)   ! Lucy potential

          R = SQRT(RSQ)
          RATIO = R/RCUT
          FORCEIJ = FORCE_FACT_LUCY*(1.0D0 - RATIO)**2

        CASE(6)   ! Soft disk potential

          FORCEIJ = 2.0d0*SOFT_N/(RSQ**(SOFT_N+1))

        CASE(7)   ! WCA fluid 
          RRR = 1.0D0/RSQ
          RR6 = RRR*RRR*RRR
          RR12 = RR6*RR6
          FORCEIJ = 48.0d0*(RR12 - 0.5d0*RR6)*RRR 

      END SELECT
   

      END SUBROUTINE
