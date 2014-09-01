      SUBROUTINE SNAP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! writes a snapshot to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double) :: SCALEFAC
      INTEGER :: I,status

      OPEN(59,FILE = 'animate.txt',STATUS = XYZ_FILE_STAT)

      SCALEFAC = 3.0D0

!!!   append frame to movie file
                                                                      
      IF(XYZ_FILE_STAT.EQ.'UNKNOWN') THEN

        readloop: DO 
          READ(59,*,IOSTAT = status)
          IF (status /= 0) EXIT
        END DO readloop 
        BACKSPACE(59)
      END IF

!!!  override replace status on file so records can be appended

      XYZ_FILE_STAT = 'UNKNOWN'

      WRITE(59,*) N
      WRITE(59,*) 'comment' 
      DO I = 1,N
            WRITE(59,1) 'N',X(I)*SCALEFAC,Y(I)*SCALEFAC,0.0D0
      END DO
      CLOSE(59)
 1    FORMAT (A2,1X,3F15.7)

      END SUBROUTINE

