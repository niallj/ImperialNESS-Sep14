      SUBROUTINE PBCS(RXX,RYY,RIX,RIY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  periodic boundary condition for non-orthoganal boxes 
!!  work in fractional space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RXX,RYY,RIX,RIY,RAI,RBI 

      CALL CARTFRAC(RXX,RYY,RAI,RBI)

      IF(RAI > 1.0d0) RAI = RAI - 1.0d0
      IF(RAI <  0.0d0) RAI = RAI + 1.0d0
      IF(RBI > 1.0d0) RBI = RBI - 1.0d0
      IF(RBI <  0.0d0) RBI = RBI + 1.0d0

      CALL FRACCART(RAI,RBI,RIX,RIY)

      END SUBROUTINE

