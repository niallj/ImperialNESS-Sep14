      SUBROUTINE INITIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  set up for 2D lattices 
!!!!  K.Travis 07/06/12 Univ. of Sheffield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE header_file 
      IMPLICIT NONE

      REAL(kind = double), DIMENSION(N) :: ATEMP,BTEMP,CTEMP 
      REAL(kind = double), DIMENSION(4) :: XTEMP,YTEMP,ZTEMP 
      REAL(kind = double) :: VX0,VY0,VZ0,DEGTORAD
      REAL(kind = double) :: RA,RB,RC,RXX,RYY,RZZ,SPACING 
      REAL(kind = double) :: XSUM,YSUM,ZSUM,PSQSUM,TSFAC,T0 
      REAL(kind = double) :: FACT,VAR,R1,R2,RANDNUM 
      INTEGER             :: IX,IY,IZ,I,J,K,NX,NY,NZ,MM,M,IA,IB,IC,INTX,INTY 

      DEGTORAD = PI/180.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! we'll set this up for an fcc lattice for now, but it should be easy
!! to change to any type of crystal provided a unit cell is read in
!! along with lattice parameters and angles. Can also use asymmetric unit 
!! if we create a new subroutine that applies the symmetry operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NX = NINT((DBLE(N/4))**(1.0d0/3.0d0))   ! atoms along each row/column

      ALEN = (DBLE(N)/DENSITY)**(1.0d0/3.0d0)
      BLEN = ALEN
      CLEN = ALEN
      ALPHA = 90.0*DEGTORAD
      BETA  = 90.0*DEGTORAD
      GAMMA = 90.0*DEGTORAD
 
      INTX = 0
      INTY = 0
      XYZ_FILE_STAT = 'REPLACE'

      CALL CELLVECT
      WRITE(*,*) 'Calculated density = ',N/(EFC(1,1)*EFC(2,2)*EFC(3,3)) 

      WRITE(*,*) 'Building fcc lattice with NX = ',NX

!!! unit cell fractional coordinates
!!! includes an offset so that atoms are not alinged along a cell boundary

      XTEMP(1) = 0.25D0
      XTEMP(2) = 0.25D0 + 0.5d0
      XTEMP(3) = 0.25D0
      XTEMP(4) = 0.25D0 + 0.5d0

      YTEMP(1) = 0.25D0
      YTEMP(2) = 0.25D0 + 0.5d0
      YTEMP(3) = 0.25D0 + 0.5d0
      YTEMP(4) = 0.25D0 

      ZTEMP(1) = 0.25D0
      ZTEMP(2) = 0.25D0
      ZTEMP(3) = 0.25D0 + 0.5d0
      ZTEMP(4) = 0.25D0 + 0.5d0 

!!!  apply translational operator to construct the supercell

      M = 0
      DO I = 1,NX
        DO J = 1,NX
          DO K = 1,NX
            DO MM = 1,4
              ATEMP(MM + M) = (XTEMP(MM) + (K - 1))/DBLE(NX)
              BTEMP(MM + M) = (YTEMP(MM) + (J - 1))/DBLE(NX)
              CTEMP(MM + M) = (ZTEMP(MM) + (I - 1))/DBLE(NX)
            END DO
            M = M + 4
          END DO
        END DO
      END DO


!!  now need to change them to cartesian coordinates

     DO I = 1,N
       RA = ATEMP(I) 
       RB = BTEMP(I) 
       RC = CTEMP(I) 
       CALL FRACCART(RA,RB,RC,RXX,RYY,RZZ)
       X(I) = RXX
       Y(I) = RYY
       Z(I) = RZZ
     END DO

!!!  assign velocities from a maxwell-boltzmann distribution
      DO I = 1,N
        MASS(I) = 1.0D00   !! assumes unit mass (can easily relax this)
        CALL MAXWELL(VX0,VY0,VZ0,TREQ,INTX,INTY)
        PX(I) = MASS(I)*VX0 
        PY(I) = MASS(I)*VY0 
        PZ(I) = MASS(I)*VZ0 
      END DO

!!!   zero the total linear momentum
      XSUM = 0.0D0
      YSUM = 0.0D0
      ZSUM = 0.0D0
      DO I = 1,N
        XSUM = XSUM + PX(I)
        YSUM = YSUM + PY(I)
        ZSUM = ZSUM + PZ(I)
      END DO
      WRITE(*,*) 'X-compt momentum at start = ',XSUM/N
      WRITE(*,*) 'Y-compt momentum at start = ',YSUM/N
      WRITE(*,*) 'Z-compt momentum at start = ',ZSUM/N
      PSQSUM = 0.0d0
      DO I = 1,N
        PX(I) = PX(I) - XSUM/N
        PY(I) = PY(I) - YSUM/N
        PZ(I) = PZ(I) - ZSUM/N
        PSQSUM = PSQSUM + (PX(I)**2 + PY(I)**2 + PZ(I)**2)/MASS(I)
      END DO
      IF(ITHERM == 0) THEN
        DEGFREE = 3.0D0*N - 3.0d0
      ELSE
        DEGFREE = 3.0D0*N - 4.0d0
      END IF
      T0 = PSQSUM/DEGFREE
      WRITE(*,*) 'Initial temperature = ',T0

!!!  rescale the temperature
      TSFAC = SQRT(TREQ/T0) 
      DO I = 1,N
        PX(I) = PX(I)*TSFAC
        PY(I) = PY(I)*TSFAC
        PZ(I) = PZ(I)*TSFAC
      END DO 

      END SUBROUTINE
