      SUBROUTINE WRITE_CONFIG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! writes out a restart file
!!!!! K.Travis 30/10/08 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      INTEGER :: I

      OPEN(93,FILE = 'restart.txt')

      DO I = 1,N
        WRITE(93,*) X(I),Y(I),Z(I)
        WRITE(93,*) PX(I),PY(I),PZ(I)
        WRITE(93,*) MASS(I)
      END DO
      WRITE(93,*) EFC(1,1), EFC(1,2), EFC(1,3)
      WRITE(93,*) EFC(2,1), EFC(2,2), EFC(2,3)
      WRITE(93,*) EFC(3,1), EFC(3,2), EFC(3,3)
      WRITE(93,*) ECF(1,1), ECF(1,2), ECF(1,3)
      WRITE(93,*) ECF(2,1), ECF(2,2), ECF(2,3)
      WRITE(93,*) ECF(3,1), ECF(3,2), ECF(3,3)
      WRITE(93,*) ALEN, BLEN, CLEN
      WRITE(93,*) ALPHA,BETA,GAMMA
      WRITE(93,*) DELTA 
	 
      CLOSE(93)

      END SUBROUTINE

