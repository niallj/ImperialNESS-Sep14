      SUBROUTINE RK4(NEQ,YY,YYP,DT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! fourth order Runge-Kutta algorithm from Hoover's SPAM book
!!!! YYN hold intermediate positions and velocities
!!!! YPN hold derivatives at intermediate times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      INTEGER :: I,NEQ
      INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)
      REAL(kind = double) :: DT
      REAL(kind = double), DIMENSION(NEQ) :: YY,YYP
      REAL(kind = double), DIMENSION(NEQ) :: YY1,YY2,YY3,YY4
      REAL(kind = double), DIMENSION(NEQ) :: YP1,YP2,YP3,YP4

      DO I = 1,NEQ
        YY1(I) = YY(I)
      END DO
      CALL RHS(YY1,YP1)

      DO I = 1,NEQ
        YY2(I) = YY(I) + 0.5D00*DT*YP1(I)
      END DO
      CALL RHS(YY2,YP2)

      DO I = 1,NEQ
        YY3(I) = YY(I) + 0.5D00*DT*YP2(I)
      END DO
      CALL RHS(YY3,YP3)

      DO I = 1,NEQ
        YY4(I) = YY(I) + DT*YP3(I)
      END DO
      CALL RHS(YY4,YP4)

      DO I = 1,NEQ
        YYP(I) = (YP1(I) + 2.0D00*(YP2(I)+YP3(I))+YP4(I))/6.0D00
      END DO

      DO I = 1,NEQ
         YY(I) = YY(I) + DT*YYP(I)
      END DO
 
      END SUBROUTINE

