      SUBROUTINE FRACCART(RA,RB,RX11,RY11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  converts fractional to cartesian cooridates 
!!!!  NB This routine automatically shifts the cartessians to have their origin at
!!!!  the cell centre
!!!!  K.Travis 06/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: RAT,RBT,RA,RB,RX11,RY11

      RAT = RA - 0.5D0
      RBT = RB - 0.5D0
      RX11 = RAT*EFC(1,1) + RBT*EFC(1,2)
      RY11 = RAT*EFC(2,1) + RBT*EFC(2,2)

      END SUBROUTINE
