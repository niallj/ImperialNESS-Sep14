      SUBROUTINE INITIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  set up for 2D lattices 
!!!!  K.Travis 05/05/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file 
      IMPLICIT NONE

      REAL(kind = double), DIMENSION(N) :: ATEMP,BTEMP 
      REAL(kind = double) :: VX0,VY0,XLEN,XLEN2,YLEN,YLEN2
      REAL(kind = double) :: RA,RB,RXX,RYY,SPACING 
      REAL(kind = double) :: XSUM,YSUM,PSQSUM,TSFAC,T0 
      REAL(kind = double) :: FACT,VAR,R1,R2,RANDNUM 
      INTEGER             :: IX,IY,I, J, NX, NY,MM,M,IA,IB,INTX,INTY 

      NX = NINT(SQRT(DBLE(N)))   ! atoms along each row/column

      IF(LATTICE_TYPE == 1) THEN
        GAMMA = PI/2.0D0
        ALEN = SQRT(DBLE(N)/DENSITY) 
        BLEN = ALEN
      ELSE IF (LATTICE_TYPE == 2) THEN
        GAMMA = PI/3.0D0
        SPACING = (1.0d0/SQRT(DENSITY))*SQRT(2.0d0)/(3.0d0**0.25d0)
        ALEN = SPACING*DBLE(NX)
        BLEN = ALEN
      END IF

      INTX = 0
      INTY = 0
      XYZ_FILE_STAT = 'REPLACE'

      CALL CELLVECT
      WRITE(*,*) 'Calculated spacing = ',SPACING 
      WRITE(*,*) 'Calculated density = ',N/(EFC(1,1)*EFC(2,2)) 

      WRITE(*,*) 'Building lattice with NX = ',NX

!! create a regular unit spaced grid using fractional coords along unit cell vectors as basis
      M = 0
      DO IA = 0,NX-1
        DO IB = 0,NX-1
           M = M + 1
           ATEMP(M)  = (DBLE(IA) + 0.5D0)/DBLE(NX) 
           BTEMP(M)  = (DBLE(IB) + 0.5D0)/DBLE(NX) 
          END DO
        END DO
      WRITE(*,*) 'Successfully constructed a regular grid!'

!!  now need to change them to cartesian coordinates

     DO I = 1,N
       RA = ATEMP(I) 
       RB = BTEMP(I) 
       CALL FRACCART(RA,RB,RXX,RYY)
       X(I) = RXX
       Y(I) = RYY
     END DO

!!!  assign velocities from a maxwell-boltzmann distribution
      DO I = 1,N
        MASS(I) = 1.0D00   !! assumes unit mass (can easily relax this)
        R1 = RANDNUM(INTX,INTY)
        R2 = RANDNUM(INTX,INTY)
        VAR = TREQ/MASS(I)
        FACT = SQRT(-2.0D0*LOG(R1)*VAR)
        VX0 = FACT*COS(2.0D0*PI*R2)
        VY0 = FACT*SIN(2.0D0*PI*R2)
        PX(I) = MASS(I)*VX0 
        PY(I) = MASS(I)*VY0 
      END DO

!!!   zero the total linear momentum
      XSUM = 0.0D0
      YSUM = 0.0D0
      DO I = 1,N
        XSUM = XSUM + PX(I)
        YSUM = YSUM + PY(I)
      END DO
      WRITE(*,*) 'X-compt momentum at start = ',XSUM/N
      WRITE(*,*) 'Y-compt momentum at start = ',YSUM/N
      PSQSUM = 0.0d0
      DO I = 1,N
        PX(I) = PX(I) - XSUM/N
        PY(I) = PY(I) - YSUM/N
        PSQSUM = PSQSUM + (PX(I)**2 + PY(I)**2)/MASS(I)
      END DO
      IF(ITHERM == 0) THEN
        DEGFREE = 2.0D0*N - 2.0d0
      ELSE
        DEGFREE = 2.0D0*N - 3.0d0
      END IF
      T0 = PSQSUM/DEGFREE
      WRITE(*,*) 'INitial temperature = ',T0

!!!  rescale the temperature
      TSFAC = SQRT(TREQ/T0) 
      DO I = 1,N
        PX(I) = PX(I)*TSFAC
        PY(I) = PY(I)*TSFAC
      END DO 

      END SUBROUTINE

      FUNCTION RANDNUM(INTX,INTY)
      IMPLICIT NONE
      INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)
      INTEGER :: I,J,INTX,INTY
      REAL(kind = double) :: RANDNUM 
      I = 1029*INTX + 1731
      J = I + 1029*INTY + 507*INTX - 1731
      INTX = MOD(I,2048)
      J = J + (I - INTX)/2048
      INTY = MOD(J,2048)
      RANDNUM = (INTX + 2048*INTY)/4194304.0D00

      END FUNCTION

