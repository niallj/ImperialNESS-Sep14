       SUBROUTINE CELLMAP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! sets up the 3D linkcell map
!!!! K.Travis 06/06/12  Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file
      IMPLICIT NONE

      INTEGER   :: IX,IY,IZ,ICELL,IMAP

!!!!! this version assumes no periodic boundaries!! 

      DO IZ = 1,NCELLZ 
        DO IY = 1,NCELLY 
          DO IX = 1,NCELLX 

            IMAP = (ICELL(IX, IY,IZ,NCELLX,NCELLY,NCELLZ) - 1)*13 

            MAP(IMAP + 1 ) = ICELL(IX + 1, IY    , IZ    , NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 2 ) = ICELL(IX + 1, IY + 1, IZ    , NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 3 ) = ICELL(IX    , IY + 1, IZ    , NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 4 ) = ICELL(IX - 1, IY + 1, IZ    , NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 5 ) = ICELL(IX + 1, IY    , IZ - 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 6 ) = ICELL(IX + 1, IY + 1, IZ - 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 7 ) = ICELL(IX    , IY + 1, IZ - 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 8 ) = ICELL(IX - 1, IY + 1, IZ - 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 9 ) = ICELL(IX + 1, IY    , IZ + 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 10) = ICELL(IX + 1, IY + 1, IZ + 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 11) = ICELL(IX    , IY + 1, IZ + 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 12) = ICELL(IX - 1, IY + 1, IZ + 1, NCELLX,NCELLY,NCELLZ)
            MAP(IMAP + 13) = ICELL(IX    , IY    , IZ + 1, NCELLX,NCELLY,NCELLZ)

          END DO
        END DO
      END DO


      END SUBROUTINE

      INTEGER FUNCTION ICELL(IXX,IYY,IZZ,NX,NY,NZ)

      IMPLICIT NONE
! declare input variables
      INTEGER,  INTENT(IN) :: IXX
      INTEGER,  INTENT(IN) :: IYY
      INTEGER,  INTENT(IN) :: IZZ
      INTEGER,  INTENT(IN) :: NX
      INTEGER,  INTENT(IN) :: NY
      INTEGER,  INTENT(IN) :: NZ
! evaluate cell index
      ICELL = 1 + MOD(IXX - 1 + NX, NX)       + &
 &                MOD(IYY - 1 + NY, NY)*NX    + &
 &                MOD(IZZ - 1 + NZ, NZ)*NX*NY 

      END FUNCTION 

