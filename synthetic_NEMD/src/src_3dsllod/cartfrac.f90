      SUBROUTINE CARTFRAC(RX11,RY11,RZ11,RA,RB,RC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  converts cartesian to fractional cooridates
!!!!  K.Travis 06/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RA,RB,RC,RX11,RY11,RZ11

      RA = RX11*ECF(1,1) + RY11*ECF(1,2) + RZ11*ECF(1,3)
      RB = RX11*ECF(2,1) + RY11*ECF(2,2) + RZ11*ECF(2,3)
      RC = RX11*ECF(3,1) + RY11*ECF(3,2) + RZ11*ECF(3,3)
      RA = RA + 0.5D0
      RB = RB + 0.5D0
      RC = RC + 0.5D0

      END SUBROUTINE
