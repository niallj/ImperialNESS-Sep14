      PROGRAM sllod_2d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Main Code for basic 2D Non-equilibrium Molecular Dynamics SLLOD method 
!!!!! currently can supoport 5 
!!!!! different pair potentials
!!!!! Uses time varying cell shape in place of Lees-Edwards boundary conditions 
!!!!! Uses link cells - best for short ranged interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Karl Travis 12/05/2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE header_file   

      IMPLICIT NONE

      INTEGER :: I,IP,J,NUMC,INN,EXIT_FLAG

      REAL(kind = double) :: XIJ,YIJ,XSQ,YSQ,RSQ,R,FACTOR
      REAL(kind = double) :: PXX_TOT,PYY_TOT,PXY_TOT,DELFAC 
      REAL(kind = double) :: ETOT,TOTAL_KE,PE,EIJ,FORCEIJ 
      REAL(kind = double) :: XXR,YYR,RXX,RYY,XSUM,YSUM,RIX,RIY 
      REAL(kind = double) :: XI,YI,XJ,YJ,RAI,RBI,RAJ,RBJ

!!    arrays that are passed to RK4
      REAL(kind = double), allocatable, DIMENSION(:) :: YY, YP  

      OPEN(96,FILE = 'energies.txt')
      WRITE(96,*) 'Time',',','Kin EN',',','POT EN',',','TOT EN',',','PRESS',',','STRESS',',','ZETA'

      OPEN(55,FILE   = 'properties',FORM   = 'unformatted',STATUS = 'unknown')

      DO I = 1,NUMPROPS
        BLOCKTOTAL(I) = 0.0d0
      END DO

!! read in the input parameters and allocate some arrays

      CALL READ_INPUTDATA  ! get the run parameters and system size(number of atoms)

      NEQ = EQN_NUMBER*N

      ALLOCATE(YY(NEQ), YP(NEQ))

      ALLOCATE(MASS(N))
      ALLOCATE(X(N),Y(N))
      ALLOCATE(PX(N),PY(N))

      IF(ISTART.EQ.0) THEN
        CALL INITIAL    ! start from a lattice
      ELSE
        CALL READ_CONFIG  ! start from a previous configuration 
      END IF
 
      IF(ITHERM == 0) THEN
        DEGFREE = 2.0D0*N - 2.0d0 ! -2 is from linear momentum conservation
      ELSE
        DEGFREE = 2.0D0*N - 3.0d0 ! extra -1 from thermostat which conserves sum psq
      END IF

!!! write out coordinates of starting configuration for debugging purposes

      OPEN(4,FILE = 'initial_config.txt')
      WRITE(4,*) 'X-posn',',','Y-posn'
      DO I = 1, N
        WRITE(4,41) X(I), Y(I)
      END DO
 41   FORMAT(F13.8,',',F13.8)
      CLOSE(4)

!!!! cell code stuff

      NCELLX  = INT(ALEN/RCUT)
      NCELLY  = INT(BLEN/RCUT)
      NCELL   = NCELLX*NCELLY
      MAPSIZE = 4*NCELL

      WRITE(*,*) 'RCUT = ',RCUT
      WRITE(*,*) 'ALEN,BLEN = ',ALEN,BLEN
      WRITE(*,*) 'Cell stats: ','ncellx = ',ncellx
      WRITE(*,*) '          : ','ncelly = ',ncelly
      WRITE(*,*) '          : ','ncell  = ',ncell
      WRITE(*,*) '          : ','mapsize  = ',mapsize

      ALLOCATE(HEAD(NCELL))
      ALLOCATE(LIST(N))
      ALLOCATE(MAP(MAPSIZE))

      NUMC = (N/NCELL) + 1
      MAXPAIR = NINT(4.5*N*NUMC)
      WRITE(*,*) 'initial estimate of MAXPAIR = ',MAXPAIR

      ALLOCATE(NI(MAXPAIR))
      ALLOCATE(NJ(MAXPAIR))

      CALL CELLMAP
      CALL SETUP_LINKCELL
      CALL CELLPAIRS(EXIT_FLAG)

      DELFAC = EFC(1,1)*SHEAR*DT

      time_loop: DO ITIME = 1, NTS

         DELTA = DELTA + DELFAC 
         IF(DELTA >= EFC(1,1)/2.) THEN
           DELTA = DELTA - EFC(1,1) 
         END IF
         CALL CELLVECT

         TIME = TIME + DT

        DO I = 1,N
          YY(I)       =  X(I)
          YY(I+N)     =  Y(I)
          YY(I+N+N)   =  PX(I)
          YY(I+N+N+N) =  PY(I)
        END DO

        CALL RK4(NEQ,YY,YP,DT)  !!! 4th order Runge Kutta integrator scheme

        DO I = 1,N
          X(I)        =  YY(I)
          Y(I)        =  YY(I+N)
          PX(I)       =  YY(I+N+N)   
          PY(I)       =  YY(I+N+N+N)
          XXR = X(I)
          YYR = Y(I)
          CALL PBCS(XXR,YYR,RXX,RYY)
          X(I) = RXX
          Y(I) = RYY
        END DO

        CALL SETUP_LINKCELL

        DO
          CALL CELLPAIRS(EXIT_FLAG)
          IF (EXIT_FLAG == 1) EXIT
          DEALLOCATE(NI)
          DEALLOCATE(NJ)
          MAXPAIR = INT(MAXPAIR*1.1) !increase maxpairs by 10% until its large enough to continue
          ALLOCATE(NI(MAXPAIR))
          ALLOCATE(NJ(MAXPAIR))
        END DO

!!!   calculation of various properties

        PE   = 0.0D0
        PXX  = 0.0D0
        PXY  = 0.0D0
        PYY  = 0.0D0
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
          XSQ = XIJ*XIJ
          YSQ = YIJ*YIJ
          RSQ = XSQ + YSQ
          IF(RSQ < CUTSQ) THEN
            CALL PAIR_ENERGY(RSQ,EIJ)
            CALL PAIR_FORCE(RSQ,FORCEIJ)
          ELSE
            EIJ = 0.0d0
            FORCEIJ = 0.0d0
          END IF
            PE =  PE  + EIJ 
            PXX = PXX + XIJ*XIJ*FORCEIJ
            PXY = PXY + XIJ*YIJ*FORCEIJ
            PYY = PYY + YIJ*YIJ*FORCEIJ
        END DO

        PE     = PE/DBLE(N)
        PV_POT = (PXX + PYY)/DBLE(2*N)

        TOTAL_KE = 0.0D0
        PKXX     = 0.0d0
        PKXY     = 0.0d0
        PKYY     = 0.0d0
        DO I = 1,N
          TOTAL_KE = TOTAL_KE + PX(I)**2/MASS(I)     + PY(I)**2/MASS(I)
          PKXX     = PKXX     + PX(I)*PX(I)/MASS(I)
          PKXY     = PKXY     + PX(I)*PY(I)/MASS(I)
          PKYY     = PKYY     + PY(I)*PY(I)/MASS(I)
        END DO
        PXX_TOT = PXX + PKXX 
        PYY_TOT = PYY + PKYY 
        PXY_TOT = PXY + PKXY 

        TKIN     = TOTAL_KE/DEGFREE                   !! ideal gas kinetic temperature
        TOTAL_KE = 0.5D0*TOTAL_KE/(DBLE(N))   !! KE per particle (= temp in 2D) 
        PV_KIN    = (PKXX + PKYY)/DBLE(2*N)

        VOLUME = EFC(1,1)*EFC(2,2)
        PRESSURE = (PKXX + PXX + PKYY + PYY)/(2.0d0*VOLUME)
        DENSITY  = DBLE(N)/VOLUME

        ETOT = TOTAL_KE + PE 

!!  pack values into property array

        PROPERTY(1) = TKIN 
        PROPERTY(2) = DENSITY 
        PROPERTY(3) = PRESSURE 
        PROPERTY(4) = TOTAL_KE 
        PROPERTY(5) = PE 
        PROPERTY(6) = ETOT 
        PROPERTY(7) = PV_KIN 
        PROPERTY(8) = PV_POT 
        PROPERTY(9) = PV_KIN + PV_POT 
        PROPERTY(10) = PXX/VOLUME 
        PROPERTY(11) = PYY/VOLUME 
        PROPERTY(12) = PXY/VOLUME 
        PROPERTY(13) = PKXX/VOLUME 
        PROPERTY(14) = PKYY/VOLUME 
        PROPERTY(15) = PKXY/VOLUME 
        PROPERTY(16) = PXX_TOT/VOLUME 
        PROPERTY(17) = PYY_TOT/VOLUME 
        PROPERTY(18) = PXY_TOT/VOLUME 
        PROPERTY(19) = SHEAR 
        PROPERTY(20) = ZETA 
        PROPERTY(21) = VOLUME 

!!  update block average accumulators

        DO I = 1,NUMPROPS
          BLOCKTOTAL(I) = BLOCKTOTAL(I) + PROPERTY(I)
        END DO

        IF(MOD(ITIME,ISUB).EQ.0) THEN

!!!  write some stuff to screen

          WRITE(*,*) '       Step ','            KE','                   PE','                     Tot Energy'
          WRITE(*,*) ITIME,TOTAL_KE,PE, ETOT

!!! write some stuff out to file as a time sequence

          WRITE(96,13) ITIME*DT,TOTAL_KE,PE,ETOT,PRESSURE,PXY_TOT/VOLUME,ZETA

 13   FORMAT(7(E18.8,','))

          CALL AVERAGER

        END IF

!!!   write xyz snapshot in concatenated file for use in animation

        IF (MOD(ITIME,MOVWRITE) == 0) THEN
          CALL SNAP
        END IF

      END DO time_loop

      OPEN(98,FILE = 'final_config.txt')
      WRITE(98,*) 'X-posn',',','Y-posn'

      DO I = 1,N
        WRITE(98,41) X(I),Y(I)
      END DO
      CLOSE(98)

!!!   compute the total linear momentum per atom at end of run

      XSUM = 0.0D0
      YSUM = 0.0D0
      DO I = 1,N
        XSUM = XSUM + PX(I)
        YSUM = YSUM + PY(I)
      END DO
      WRITE(*,*) 'X-compt momentum at end = ',XSUM/N
      WRITE(*,*) 'Y-compt momentum at end = ',YSUM/N

      CALL ANALYSE_PROPERTIES



      CALL WRITE_CONFIG

      CLOSE(96)


      END PROGRAM 

