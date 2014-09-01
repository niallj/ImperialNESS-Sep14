       SUBROUTINE CELLMAP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! sets up the 2D linkcell map
!!!! K.Travis 24/02/12  Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      INTEGER   :: IX,IY,ICELL,IMAP

      DO IY = 1,NCELLY 
        DO IX = 1,NCELLX 

          IMAP = (ICELL(IX, IY,NCELLX,NCELLY) - 1)*4 

          MAP(IMAP + 1) = ICELL(IX + 1, IY    , NCELLX, NCELLY)
          MAP(IMAP + 2) = ICELL(IX + 1, IY + 1, NCELLX, NCELLY)
          MAP(IMAP + 3) = ICELL(IX    , IY + 1, NCELLX, NCELLY)
          MAP(IMAP + 4) = ICELL(IX - 1, IY + 1, NCELLX, NCELLY)

        END DO
      END DO


      END SUBROUTINE

      INTEGER FUNCTION ICELL(IXX,IYY,NX,NY)
      IMPLICIT NONE
! declare input variables
      INTEGER,  INTENT(IN) :: IXX
      INTEGER,  INTENT(IN) :: IYY
      INTEGER,  INTENT(IN) :: NX
      INTEGER,  INTENT(IN) :: NY

! evaluate cell index

      ICELL = 1 + MOD(IXX - 1 + NX,NX) + MOD(IYY - 1 + NY, NY)*NX 

      END FUNCTION 

