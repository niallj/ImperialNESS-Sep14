SUBROUTINE READ_CONFIG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Reads in a restart file 
!!!!  K.Travis 24/03/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE header_file
IMPLICIT NONE
     
REAL(KIND = DOUBLE) :: DUMMY,DET
INTEGER :: I

OPEN(4,FILE = 'restart.txt')

DO I = 1,N
   READ(4,*) X(I),Y(I)
   READ(4,*) PX(I),PY(I)
   READ(4,*) MASS(I) 
END DO
READ(4,*) EFC(1,1), EFC(1,2) 
READ(4,*) EFC(2,1), EFC(2,2) 
READ(4,*) ECF(1,1), ECF(1,2)
READ(4,*) ECF(2,1), ECF(2,2)
READ(4,*) ALEN,BLEN,GAMMA, DELTA 

VOLUME = EFC(1,1)*EFC(2,2) 

XYZ_FILE_STAT = 'UNKNOWN'

CLOSE(4)

END SUBROUTINE

