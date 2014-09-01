      SUBROUTINE PBCS(RXX,RYY,RZZ,RIX,RIY,RIZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  periodic boundary condition for non-orthoganal boxes 
!!  work in fractional space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RXX,RYY,RZZ,RIX,RIY,RIZ,RAI,RBI,RCI 

      CALL CARTFRAC(RXX,RYY,RZZ,RAI,RBI,RCI)

      IF(RAI > 1.0d0) RAI = RAI - 1.0d0
      IF(RAI <  0.0d0) RAI = RAI + 1.0d0
      IF(RBI > 1.0d0) RBI = RBI - 1.0d0
      IF(RBI <  0.0d0) RBI = RBI + 1.0d0
      IF(RCI > 1.0d0) RCI = RCI - 1.0d0
      IF(RCI <  0.0d0) RCI = RCI + 1.0d0

      CALL FRACCART(RAI,RBI,RCI,RIX,RIY,RIZ)

      END SUBROUTINE

