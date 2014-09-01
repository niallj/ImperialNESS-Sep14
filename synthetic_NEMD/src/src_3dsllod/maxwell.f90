        SUBROUTINE MAXWELL(VX,VY,VZ,TEMP,INTX,INTY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  Subroutine to select two cartesian components of     !!!!! 
!!!!!  velocity from a Maxwell-Boltzmann distribution at    !!!!! 
!!!!!  temperature, TEMP_REQ.                               !!!!! 
!!!!!  uses the Box-Muller method.                          !!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)
        REAL(kind = double), PARAMETER :: PI2 = 2.0d0*3.141592654D0
        REAL(kind = double) :: VAR,FACT,FACT2,VX,VY,VZ,TEMP
        REAL(kind = double) :: R1,R2,R3,R4,RANDNUM 
        INTEGER :: INTX,INTY

        R1 = RANDNUM(INTX,INTY)
        R2 = RANDNUM(INTX,INTY)
        R3 = RANDNUM(INTX,INTY)
        R4 = RANDNUM(INTX,INTY)

        VAR = TEMP

        FACT = SQRT(-2.0d0*LOG(R1)*VAR)
        VX = FACT*COS(PI2*R2)
        VY = FACT*SIN(PI2*R2)

        FACT2 = SQRT(-2.0d0*LOG(R3)*VAR)
        VZ = FACT2*COS(PI2*R4)

        END SUBROUTINE

      FUNCTION RANDNUM(INTX,INTY)

      IMPLICIT NONE
      INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)
      INTEGER :: I,J,INTX,INTY
      REAL(kind = double) :: RANDNUM
      I = 1029*INTX + 1731
      J = I + 1029*INTY + 507*INTX - 1731
      INTX = MOD(I,2048)
      J = J + (I - INTX)/2048
      INTY = MOD(J,2048)
      RANDNUM = (INTX + 2048*INTY)/4194304.0D00

      END FUNCTION



