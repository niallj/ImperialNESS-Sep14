      SUBROUTINE RHS(YY,YP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    01/09/14, K.Travis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file 
      IMPLICIT NONE

      INTEGER :: I, J, IP 
      REAL(kind = double) :: XI,XJ,ZJ,YI,YJ,ZI,XIJ,YIJ,ZIJ
      REAL(kind = double) :: R, RSQ, FACTOR, FORCE_MN,DELR
      REAL(kind = double) :: RATIO, FORCE_REP,SUMPXPY,TFLAG
      REAL(kind = double) :: RRR,RR6,RR12,FORCEIJ,PXYP,PXYK,PXYV 
      REAL(kind = double) :: RDIF,RDIF2,RDIF3,RXX,RYY,RZZ,RAI,RBI,RCI,RAJ,RBJ,RCJ 
      REAL(kind = double) :: ESUM,TK,SUMFP,SUMPP,SUMZFX 

      REAL(kind = double), DIMENSION(NEQ) ::  YY,YP
      REAL(kind = double), DIMENSION(N)   ::  FX,FY,FZ

      FX = (/ (I*0.0D0,I = 1,N) /)
      FY = (/ (I*0.0D0,I = 1,N) /)
      FZ = (/ (I*0.0D0,I = 1,N) /)

!     store arguments into local arrays

      DO I = 1,N
        X(I)    = YY(I)
        Y(I)    = YY(I+N)
        Z(I)    = YY(I+N+N)
        PX(I)   = YY(I+N+N+N)
        PY(I)   = YY(I+N+N+N+N)
        PZ(I)   = YY(I+N+N+N+N+N)
      END DO
      ZETA      = YY(N+N+N+N+N+N+1)


!     compute the instantaneous temperature

      ESUM = 0.0D0
      DO I = 1,N
          ESUM = ESUM + PX(I)**2/MASS(I) + PY(I)**2/MASS(I) + PZ(I)**2/MASS(I)
      END DO

      TK = ESUM/DEGFREE

! calculate pair forces
      DO IP = 1, NPAIRS
        I = NI(IP)
        J = NJ(IP)
        XI = X(I)
        YI = Y(I)
        ZI = Z(I)
        CALL CARTFRAC(XI,YI,ZI,RAI,RBI,RCI)
        XJ = X(J)
        YJ = Y(J)
        ZJ = Z(J)
        CALL CARTFRAC(XJ,YJ,ZJ,RAJ,RBJ,RCJ)
        RAJ = RAJ - ANINT(RAJ-RAI)
        RBJ = RBJ - ANINT(RBJ-RBI)
        RCJ = RCJ - ANINT(RCJ-RCI)
        CALL FRACCART(RAJ,RBJ,RCJ,XJ,YJ,ZJ)
        XIJ = XI - XJ
        YIJ = YI - YJ
        ZIJ = ZI - ZJ
        RSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ 
        IF(RSQ < CUTSQ) THEN

          CALL PAIR_FORCE(RSQ,FORCEIJ) 
          FX(I) = FX(I) + FORCEIJ*XIJ
          FY(I) = FY(I) + FORCEIJ*YIJ
          FZ(I) = FZ(I) + FORCEIJ*ZIJ
          FX(J) = FX(J) - FORCEIJ*XIJ
          FY(J) = FY(J) - FORCEIJ*YIJ
          FZ(J) = FZ(J) - FORCEIJ*ZIJ
       END IF

      END DO
      
      IF(ITHERM.EQ.1) THEN
        TFLAG = 1.0D0
      ELSE
        TFLAG = 0.0D0
      END IF

 !!! pack all first derivatives into a 1-D array

      DO I = 1, N
          YP(I)                  = PX(I)/MASS(I) + SHEAR*Y(I) 
          YP(I+N)                = PY(I)/MASS(I) 
          YP(I+N+N)              = PZ(I)/MASS(I) 
          YP(I+N+N+N)            = FX(I) - ZETA*PX(I) - SHEAR*PY(I) 
          YP(I+N+N+N+N)          = FY(I) - ZETA*PY(I) 
          YP(I+N+N+N+N+N)        = FZ(I) - ZETA*PZ(I) 
      END DO
      YP(N+N+N+N+N+N+1)          = TFLAG*(TK - TREQ)/(TAU_TEMP_SQD*TREQ)

      END SUBROUTINE

