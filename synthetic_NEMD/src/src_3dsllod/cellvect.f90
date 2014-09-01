      SUBROUTINE CELLVECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generates the vectors used to convert fractionals to
!! Cartesians and cartesians to fractionals
!! Using the cell parameters a,b 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE
 
      REAL(kind = double) :: DET
      REAL(kind = double), DIMENSION(3,3) :: DUMMY1,DUMMY2 
      INTEGER :: I,J

!!! EA is on x-axis so is (a,0,0)

      EFC(1,1) = ALEN
      EFC(1,2) = BLEN*COS(GAMMA) + DELTA 
      EFC(1,3) = CLEN*COS(GAMMA) 

!!! EB lies in xy plane at an angle gamma to x

      EFC(2,1) = 0.0D0
      EFC(2,2) = BLEN*SIN(GAMMA) 
      EFC(2,3) = CLEN*(COS(ALPHA) -(COS(GAMMA)*COS(BETA)))/SIN(GAMMA) 

!!! EC points into space and has components in all directions

      EFC(3,1) = 0.0D0
      EFC(3,2) = 0.0D0 
      EFC(3,3) = CLEN - CLEN*COS(BETA)*COS(BETA) +                        & 
      &           CLEN*((COS(ALPHA) - COS(ALPHA)*COS(BETA)) / SIN(GAMMA)) &
      &               *((COS(ALPHA)-COS(ALPHA)*COS(BETA))/SIN(GAMMA)) 

!!! copy EFC array elements to dummy1 matrix

        DO I = 1,3
          DO J = 1,3
            DUMMY1(I,J) = EFC(I,J)
          END DO
        END DO

        CALL INVERT(DUMMY1,DUMMY2)

!!! copy dummy2 array elements to ECF matrix

        DO I = 1,3
          DO J = 1,3
            ECF(I,J) = DUMMY2(I,J)
          END DO
        END DO


      END SUBROUTINE

