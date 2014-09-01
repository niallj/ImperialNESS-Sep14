      SUBROUTINE SETUP_LINKCELL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  routine to determine which link cell each particle resides in
!!!  and to form appropriate lists. This is used for an order
!!!  N force/energy evaluation 
!!!  Written by K.Travis 24/02/12 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file  !!use stuff in main module

      IMPLICIT NONE

      INTEGER :: I,ICELL,ICX,ICY

      REAL(kind = double) :: RA,RB,ACELL,BCELL,XIS,YIS,XX,YY


!!!  zero arrays

      HEAD = (/ (I*0,I = 1,NCELL) /)
      LIST = (/ (I*0,I = 1,N)     /)

!!! sort particles into cells

      DO I = 1,N
        XX = X(I)
        YY = Y(I)
        CALL CARTFRAC(XX,YY,RA,RB)
        ICX = INT(RA*NCELLX)
        ICX = MAX(0,ICX)
        ICX = MIN(ICX,NCELLX-1)
        ICY = INT(RB*NCELLY)
        ICY = MAX(0,ICY)
        ICY = MIN(ICY,NCELLY-1)
        ICELL = 1 + ICX + ICY*NCELLX
        LIST(I) = HEAD(ICELL)
        HEAD(ICELL) = I
      END DO

      END SUBROUTINE
