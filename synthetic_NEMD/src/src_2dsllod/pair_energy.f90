      SUBROUTINE PAIR_ENERGY(RSQ,EIJ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculates the energy between a pair of atoms 
!! Current supports 6 different pair potentials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: EIJ,RSQ 
      REAL(kind = double) :: RRR,RR6,RR12 
      REAL(kind = double) :: R,RDIF,RDIF2,RDIF3 
      REAL(kind = double) :: RATIO,FACT1
      REAL(kind = double) :: FACTOR

      SELECT CASE (IPOT_TYPE) ! choose option depending on potential type

        CASE(1)   ! Cut LJ potential

          RRR = 1.0D0/RSQ
          RR6 = RRR*RRR*RRR
          RR12 = RR6*RR6
          EIJ = 4.0d0*EPS*(RR12 - RR6) 

        CASE(2)   ! Splined-LJ potential

          IF(RSQ.LE.RCUT_INNER_SQ) THEN
            RRR = 1.0D0/RSQ
            RR6 = RRR*RRR*RRR
            RR12 = RR6*RR6
            EIJ = 4.0d0*EPS*(RR12 - RR6)
          ELSE
            R = SQRT(RSQ)
            RDIF = R - RCUT
            RDIF2 = RDIF*RDIF
            RDIF3 = RDIF2*RDIF
            EIJ = SPLINE_R2*RDIF2 + SPLINE_R3*RDIF3
          END IF

        CASE(3)   ! Smooth repulsive potential

          RATIO = RSQ/CUTSQ
          EIJ = EPS*(1.0D0 - RATIO)**4

        CASE(4)   ! m-n family of potentials

          FACTOR = 2.0D0 - RSQ
          EIJ = (M_RATIO*FACTOR**N_POWER - N_RATIO*FACTOR**M_POWER)

        CASE(5)   ! Lucy potential

          R = SQRT(RSQ)
          RATIO = R/RCUT
          FACT1 = (1.0D0 - RATIO)
          EIJ = EFACT_LUCY*FACT1**3*(1.0D0 + 3.0D0*RATIO)

        CASE(6)   ! Soft disk potential

          EIJ = 1.0D0/(RSQ**SOFT_N)

        CASE(7)   ! WCA potential

          RRR = 1.0D0/RSQ
          RR6 = RRR*RRR*RRR
          RR12 = RR6*RR6
          EIJ = 4.0d0*(RR12 - RR6) + 1.0d0 

      END SELECT
   

      END SUBROUTINE
