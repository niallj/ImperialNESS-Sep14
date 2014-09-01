      SUBROUTINE CELLVECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generates the vectors used to convert fractionals to
!! Cartesians and cartesians to fractionals
!! Using the cell parameters a,b 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE
 
      REAL(kind = double) :: DET

      EFC(1,1) = ALEN
      EFC(1,2) = BLEN*COS(GAMMA) + DELTA 
      EFC(2,1) = 0.0D0
      EFC(2,2) = BLEN*SIN(GAMMA) 

      DET = ABS(EFC(1,1)*EFC(2,2) - EFC(1,2)*EFC(2,1))

      ECF(1,1) = EFC(2,2)/DET
      ECF(1,2) = -EFC(1,2)/DET
      ECF(2,1) = -EFC(2,1)/DET
      ECF(2,2) = EFC(1,1)/DET

      END SUBROUTINE

