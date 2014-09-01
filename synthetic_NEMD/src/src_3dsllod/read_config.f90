SUBROUTINE READ_CONFIG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Reads in a restart file in which the fluid is a 
!!!!  circular droplet
!!!!  K.Travis 24/03/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE header_file
IMPLICIT NONE
     
REAL(KIND = DOUBLE) :: DUMMY,DET
INTEGER :: I

OPEN(4,FILE = 'restart.txt')

DO I = 1,N
   READ(4,*) X(I),Y(I),Z(I)
   READ(4,*) PX(I),PY(I),PZ(I)
   READ(4,*) MASS(I) 
END DO
READ(4,*) EFC(1,1), EFC(1,2),EFC(1,3) 
READ(4,*) EFC(2,1), EFC(2,2),EFC(2,3) 
READ(4,*) EFC(3,1), EFC(3,2),EFC(3,3) 
READ(4,*) ECF(1,1), ECF(1,2),ECF(1,3)
READ(4,*) ECF(2,1), ECF(2,2),ECF(2,3)
READ(4,*) ECF(3,1), ECF(3,2),ECF(3,3)
READ(4,*) ALEN,BLEN,CLEN 
READ(4,*) ALPHA,BETA,GAMMA
READ(4,*) DELTA 

VOLUME = EFC(1,1)*EFC(2,2)*EFC(3,3) 

XYZ_FILE_STAT = 'UNKNOWN'

CLOSE(4)

END SUBROUTINE

