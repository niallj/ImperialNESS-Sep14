      SUBROUTINE AVERAGER      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  sets up block averages 
!!!!  K.Travis 10/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      INTEGER  :: I

      DO I = 1,NUMPROPS
        BLOCKAVG(I) = BLOCKTOTAL(I)/DBLE(ISUB)
      END DO

!!  zero all block averaged accumulators

      DO I = 1,NUMPROPS
        BLOCKTOTAL(I) = 0.0D0
      END DO

!! write the averages to file in binary format (saves disk space)

      WRITE(55) ITIME,BLOCKAVG

      END SUBROUTINE


