      SUBROUTINE RHS(YY,YP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Provides right hand sides of eqns of motion
!!!   K.Travis 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file 
      IMPLICIT NONE

      INTEGER :: I, J, IP 
      REAL(kind = double) :: XI,XJ,YI,YJ,XIJ, YIJ
      REAL(kind = double) :: R, RSQ, FACTOR, FORCE_MN,DELR
      REAL(kind = double) :: RATIO, FORCE_REP,SUMPXPY
      REAL(kind = double) :: RRR,RR6,RR12,FORCEIJ,PXYP,PXYK,PXYV 
      REAL(kind = double) :: RDIF,RDIF2,RDIF3,RXX,RYY,RAI,RBI,RAJ,RBJ 
      REAL(kind = double) :: ESUM,TK,DRIFT, SUMFP,SUMPP

      REAL(kind = double), DIMENSION(NEQ) ::  YY,YP
      REAL(kind = double), DIMENSION(N)   ::  FX,FY

      FX = (/ (I*0.0D0,I = 1,N) /)
      FY = (/ (I*0.0D0,I = 1,N) /)

!     store arguments into local arrays

      DO I = 1,N
        X(I)    = YY(I)
        Y(I)    = YY(I+N)
        PX(I)   = YY(I+N+N)
        PY(I)   = YY(I+N+N+N)
      END DO


!     compute the instantaneous temperature

      ESUM = 0.0D0
      PXYK = 0.0d0
      DO I = 1,N
          ESUM = ESUM + PX(I)**2/MASS(I) + PY(I)**2/MASS(I)
          PXYK = PXYK + PX(I)*PY(I)/MASS(I)
      END DO

      TK = ESUM/DEGFREE

! calculate pair forces
      PXYP = 0.0D0
      DO IP = 1, NPAIRS
        I = NI(IP)
        J = NJ(IP)
        XI = X(I)
        YI = Y(I)
        CALL CARTFRAC(XI,YI,RAI,RBI)
        XJ = X(J)
        YJ = Y(J)
        CALL CARTFRAC(XJ,YJ,RAJ,RBJ)
        RAJ = RAJ - ANINT(RAJ-RAI)
        RBJ = RBJ - ANINT(RBJ-RBI)
        CALL FRACCART(RAJ,RBJ,XJ,YJ)
        XIJ = XI - XJ
        YIJ = YI - YJ
        RSQ = XIJ*XIJ + YIJ*YIJ 

        IF(RSQ < CUTSQ) THEN

          CALL PAIR_FORCE(RSQ,FORCEIJ) 

          FX(I) = FX(I) + FORCEIJ*XIJ
          FY(I) = FY(I) + FORCEIJ*YIJ
          FX(J) = FX(J) - FORCEIJ*XIJ
          FY(J) = FY(J) - FORCEIJ*YIJ
          PXYP = PXYP + XIJ*YIJ*FORCEIJ
       END IF

      END DO
      
      PXYV = 0.5*(PXYP + PXYK)

!! compute the Gaussian constraint thermostat multiplier

      SELECT CASE (ITHERM) ! choose option depending on potential type

      CASE(1)   ! Gaussian thermostat 

        SUMFP = 0.0D0
        SUMPP = 0.0D0
        SUMPXPY = 0.0D0
        DO I = 1,N
          SUMFP = SUMFP + FX(I)*PX(I) + FY(I)*PY(I)
          SUMPP = SUMPP + PX(I)*PX(I) + PY(I)*PY(I)
          SUMPXPY = SUMPXPY + PX(I)*PY(I)
        END DO
        ZETA = (SUMFP - SHEAR*SUMPXPY)/SUMPP
        DRIFT = BETA*(TK - TREQ)

      CASE(2)   ! Gaussian ergostat

        SUMPP = 0.0D0
        DO I = 1,N
          SUMPP = SUMPP + PX(I)*PX(I)/MASS(I) + PY(I)*PY(I)/MASS(I)
        END DO
        ZETA = -SHEAR*PXYV/SUMPP
        DRIFT = 0.0d0

      CASE DEFAULT

        ZETA = 0.0D0
        DRIFT = 0.0d0 

      END SELECT 

 !!! pack all first derivatives into a 1-D array

      DO I = 1, N
          YP(I)                  = PX(I)/MASS(I) + SHEAR*Y(I) 
          YP(I+N)                = PY(I)/MASS(I) 
          YP(I+N+N)              = FX(I) - (ZETA + DRIFT)*PX(I) - SHEAR*PY(I) 
          YP(I+N+N+N)            = FY(I) - (ZETA + DRIFT)*PY(I) 
      END DO

      END SUBROUTINE

