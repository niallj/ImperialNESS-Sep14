      SUBROUTINE ANALYSE_PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! reads in block averages from disk
!!!! computes mean, standard deviation, standard error and
!!!! statistical inefficiency time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file
      IMPLICIT NONE

      REAL(kind = double), ALLOCATABLE, DIMENSION(:,:) :: A
      REAL(kind = double), DIMENSION(50,NUMPROPS) :: SIS 
      REAL(kind = double), DIMENSION(NUMPROPS) :: SE,B, C, D, MEAN, FACT 
      REAL(kind = double), DIMENSION(50) :: XST 
      REAL(kind = double) :: RNT,S,AV,T,RNS,AVN,RMSD,CI,SI,YMX,SIMX
      INTEGER :: STATUS
      INTEGER :: I,J,K,NSUB,INT,NIT
      INTEGER :: KK,L,NMAX,NTOTAL,NTOT1 
      INTEGER :: NTIM = 0 

      CHARACTER(20), DIMENSION(NUMPROPS) :: PROPNAMES


      PROPNAMES(1) = 'kinetic Temp'
      PROPNAMES(2) = 'density'
      PROPNAMES(3) = 'pressure'
      PROPNAMES(4) = 'kinetic energy/N'
      PROPNAMES(5) = 'potential energy/N'
      PROPNAMES(6) = 'total energy/N'
      PROPNAMES(7) = 'PV/N: Kinetic'
      PROPNAMES(8) = 'PV/N: Potential'
      PROPNAMES(9) = 'PV/N: total'
      PROPNAMES(10) = 'PXX: potential '
      PROPNAMES(11) = 'PYY: potential '
      PROPNAMES(12) = 'PXY: potential '
      PROPNAMES(13) = 'PYZ: potential '
      PROPNAMES(14) = 'PZZ: potential '

      PROPNAMES(15) = 'PXX: kinetic '
      PROPNAMES(16) = 'PYY: kinetic '
      PROPNAMES(17) = 'PXY: kinetic '
      PROPNAMES(18) = 'PYZ: kinetic '
      PROPNAMES(19) = 'PZZ: kinetic '
      PROPNAMES(20) = 'PXX_TOT '
      PROPNAMES(21) = 'PYY_TOT '
      PROPNAMES(22) = 'PXY_TOT '
      PROPNAMES(23) = 'PYZ_TOT '
      PROPNAMES(24) = 'PZZ_TOT '
      PROPNAMES(25) = 'SHEAR '
      PROPNAMES(26) = 'ZETA '
      PROPNAMES(27) = 'VOLUME '
      
      OPEN(12,FILE = 'results.txt',FORM = 'formatted',STATUS = 'unknown')


!!!!  read binary file of data cccccccccccccccc

      REWIND(55)
      
      readloop: DO     
        READ(55,IOSTAT = status) ITIME,BLOCKAVG
        IF(status /= 0) EXIT
        NTIM = NTIM + 1
      END DO readloop

      readif: IF(status > 0) THEN ! a read error occured. 
         WRITE(*,10) (NTIM + 1)
         10 FORMAT ('0', 'An error occurred reading line ', I7) 
         STOP
      ELSE ! got to end of the file OK
         WRITE(*,20) NTIM
         20 FORMAT ('0', 'End of file reached. There were ', I7,&
                         ' values in the file.') 
      END IF readif

      ALLOCATE(A(NTIM,NUMPROPS))

      REWIND(55)

      DO J = 1, NTIM
        READ(55) ITIME,BLOCKAVG
        DO I = 1, NUMPROPS
          A(J,I) = BLOCKAVG(I)
        END DO 
      END DO 

      RNT = NTIM 

!!!  get averages

      DO I = 1, NUMPROPS
        S = 0.0d0
        DO J = 1, NTIM
          S = S + A(J,I)
        END DO 
        B(I) = S/RNT
      END DO 
                                                           
!!! get variances

      DO I = 1, NUMPROPS
        AV = B(I)
        S = 0.0d0
        DO J = 1, NTIM
          S = S + (A(J,I) - AV)**2
        END DO 
        C(I) = S/RNT
      END DO 

      NSUB = 10
      INT = 0
      NIT = 0
      YMX = 0.0

!!! following code written donkeys years ago by David Brown in fortran 77.
!!! its not very structured so needs an overall at some stage

  500 CONTINUE 
      NSUB = NSUB + 10*(1 + (NSUB/100))
      RNS = NSUB 
      RNS = 1.0/RNS
      NMAX = NTIM/NSUB
      NIT = NIT + 1
      IF(NMAX < 5)  GOTO 501
      IF(NIT  > 50) GOTO 501
      INT = INT + 1
      XST(INT) = NSUB 
      NTOTAL = NMAX*NSUB
      NTOT1 = NTOTAL + 1
      RNT = NTOTAL 
      DO 502 J = 1,NUMPROPS
        AVN = 0.0d0
        IF(NTOTAL == NTIM) GOTO 503
        DO 504 I = NTOT1,NTIM
          AVN = AVN + A(I,J)
  504  CONTINUE 
  503 CONTINUE 
      MEAN(J) = (NTIM*B(J) - AVN)/RNT
  502 CONTINUE 
      DO 509 J = 1,NUMPROPS
        I = 0
        T = 0.0d0
        DO 510 KK = 1,NMAX
          S = 0.0d0
          DO 511 l = 1,NSUB
            I = I + 1
            S = S + A(I,J)
  511     CONTINUE 
        S = S*RNS
        T = T + (S - MEAN(J))**2
  510 CONTINUE 
      D(J) = T/DBLE(NMAX)
  509 CONTINUE 
                                                                     
      DO I = 1,NUMPROPS
        RMSD = SQRT(D(I))
        CI = C(I)
        IF(C(I) == 0.0d0) CI = 1.0d0
        SI = NSUB*D(I)/CI
        SIS(INT,I) = SI 
        YMX = MAX(YMX,SI)
!        WRITE(12,601) PROPNAMES(I),MEAN(I),D(I),RMSD,SI
  601 format(a20,3x,4(e13.6,2x))
      END DO 
                                                                    
      GOTO 500

  501 CONTINUE 

      WRITE(12,698)
  698 format(///,'property',14x,'average',5x,'stand.error',&
     & 7x,'rmsd',8x,'stat.ineff.')
      WRITE(12,699)
  699 format('--------',14x,'-------',5x,'-----------',&
     & 7x,'----',8x,'-----------')
      WRITE(12,*)

      DO I = 1,NUMPROPS
        SIMX = 0.0d0
        DO J = 1,INT
          SIMX = MAX(SIMX,SIS(J,I))
        END DO 
        SE(I) = SQRT(SIMX*C(I)/NTIM)
        RMSD = SQRT(C(I))
        WRITE(12,702) PROPNAMES(I),B(I),SE(I),RMSD,SIMX
      END DO 
 702  FORMAT(a20,2x,4(e13.6,2x))

!!! write out some data in new file for use in perl scripts

      OPEN(22,FILE = 'Tabulated_results.txt')

      WRITE(22,*) B(25),',',SE(25)  !SHEAR RATE
      WRITE(22,*) B(1),',',SE(1)  !temperature
      WRITE(22,*) B(2),',',SE(2)  !density
      WRITE(22,*) B(27),',',SE(27)  !VOLUME
      WRITE(22,*) B(3),',',SE(3)  !pressure
      WRITE(22,*) B(6),',',SE(6)  !ETOT per particle
      WRITE(22,*) B(9),',',SE(9)  !PV per particle
      WRITE(22,*) B(20),',',SE(20)  !PXX
      WRITE(22,*) B(21),',',SE(21)  !PYY
      WRITE(22,*) B(24),',',SE(24)  !PZZ
      WRITE(22,*) B(22),',',SE(22)  !PXY
      WRITE(22,*) B(26),',',SE(26)  !ZETA
      CLOSE(22)

      END SUBROUTINE
