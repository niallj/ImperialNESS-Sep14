      SUBROUTINE CELLPAIRS(IFG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Returns a list of npair ij pairs 
!!!!! Utilises cell code for efficient sorting
!!!!! K.Travis 24/02/12 Univ. of Sheffield.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file    !! use stuff in main module named header_file

      INTEGER :: M, ICELL, I, J, NABOR, JCELLO, JCELL,IFG

      NPAIRS = 0
      IFG = 1
      DO M = 1,MAXPAIR
        NI(M) = 0
        NJ(M) = 0
      END DO

!!!  loop over all cells

      DO ICELL = 1,NCELL 
         I = HEAD(ICELL)
        DO WHILE (I.GT.0)
          J = LIST(I)
          DO WHILE (J.GT.0)
            NPAIRS = NPAIRS + 1
            IF(NPAIRS > MAXPAIR) THEN
              IFG = 0
              RETURN
            END IF
            NI(NPAIRS) = I
            NJ(NPAIRS) = J
            J = LIST(J)
          END DO

!!! loop over neighbouring cells

          JCELLO = 13*(ICELL - 1)
          DO NABOR = 1,13
            JCELL = MAP(JCELLO + NABOR)
            IF(JCELL.GT.0) THEN
              J = HEAD(JCELL)
              DO WHILE(J.NE.0)
                NPAIRS = NPAIRS + 1
                IF(NPAIRS > MAXPAIR) THEN
                  IFG = 0
                  RETURN
                END IF
                NI(NPAIRS) = I
                NJ(NPAIRS) = J
                J = LIST(J)
              END DO
            END IF
          END DO
          I = LIST(I)
        END DO
      END DO

      END SUBROUTINE

