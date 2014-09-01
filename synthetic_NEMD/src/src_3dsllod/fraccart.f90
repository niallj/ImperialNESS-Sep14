      SUBROUTINE FRACCART(RA,RB,RC,RX11,RY11,RZ11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  converts fractional to cartesian cooridates 
!!!!  NB This routine automatically shifts the cartessians to have their origin at
!!!!  the cell centre
!!!!  K.Travis 06/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RAT,RBT,RCT,RA,RB,RC,RX11,RY11,RZ11

      RAT = RA - 0.5D0
      RBT = RB - 0.5D0
      RCT = RC - 0.5D0
      RX11 = RAT*EFC(1,1) + RBT*EFC(1,2) + RCT*EFC(1,3)
      RY11 = RAT*EFC(2,1) + RBT*EFC(2,2) + RCT*EFC(2,3)
      RZ11 = RAT*EFC(3,1) + RBT*EFC(3,2) + RCT*EFC(3,3)

      END SUBROUTINE
