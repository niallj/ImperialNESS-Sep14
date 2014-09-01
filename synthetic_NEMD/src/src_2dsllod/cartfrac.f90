      SUBROUTINE CARTFRAC(RX11,RY11,RA,RB)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  converts cartesian to fractional cooridates
!!!!  K.Travis 06/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RA,RB,RX11,RY11

      RA = RX11*ECF(1,1) + RY11*ECF(1,2)
      RB = RX11*ECF(2,1) + RY11*ECF(2,2)
      RA = RA + 0.5D0
      RB = RB + 0.5D0

      END SUBROUTINE
