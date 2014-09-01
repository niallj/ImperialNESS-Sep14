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
        WRITE(93,*) X(I),Y(I)
        WRITE(93,*) PX(I),PY(I)
        WRITE(93,*) MASS(I)
      END DO
      WRITE(93,*) EFC(1,1), EFC(1,2)
      WRITE(93,*) EFC(2,1), EFC(2,2)
      WRITE(93,*) ECF(1,1), ECF(1,2)
      WRITE(93,*) ECF(2,1), ECF(2,2)
      WRITE(93,*) ALEN, BLEN, GAMMA, DELTA 
	 
      CLOSE(93)

      END SUBROUTINE

