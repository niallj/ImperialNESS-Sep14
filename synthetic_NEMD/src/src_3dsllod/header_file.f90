MODULE header_file
!
! replacement for the old fortran77 include file with common blocks
! here we declare some shared memory variables such as position and momenta arrays
! they are now allocatable, with the size being determined through the input file

IMPLICIT NONE
PUBLIC
SAVE

INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15,99)

INTEGER, PARAMETER :: EQN_NUMBER = 7 ! number of right hand side equations to be solved 
INTEGER, PARAMETER :: NUMPROPS   = 27 ! number of properties used in block averaging scheme 

REAL(kind = double), PARAMETER :: PI = 3.141592654D0
REAL(kind = double), PARAMETER :: PI2 = 2.0D0*PI 
REAL(kind = double), PARAMETER :: CHI = 10.0D0 

REAL(kind = double), DIMENSION(3,3) :: EFC, ECF 
REAL(kind = double), DIMENSION(NUMPROPS) :: BLOCKTOTAL,BLOCKAVG,PROPERTY 

REAL(kind = double), ALLOCATABLE, DIMENSION(:) :: X, Y, Z
REAL(kind = double), ALLOCATABLE, DIMENSION(:) :: PX, PY, PZ
REAL(kind = double), ALLOCATABLE, DIMENSION(:) :: MASS 

INTEGER, ALLOCATABLE, DIMENSION(:) :: NI, NJ, MAP  
INTEGER, ALLOCATABLE, DIMENSION(:) :: HEAD, LIST  

REAL(kind = double) :: ZETA,DENSITY, SIGMA, EPS, DT, VOLUME,TIME
REAL(kind = double) :: SIGSQ, EPSIG8, DIFF_MN, M_RATIO, N_RATIO 
REAL(kind = double) :: FORCE_MN_FACT, RCUT,SHEAR,DELTA,ALPHA,BETA,GAMMA 
REAL(kind = double) :: RCUT_INNER,SPLINE_R2,SPLINE_R3,RCUT_INNER_SQ
REAL(kind = double) :: CUTSQ,TREQ,ALEN,BLEN,CLEN,TIME_STRAIN
REAL(kind = double) :: SOFT_N,H_LUCY,EFACT_LUCY,FORCE_FACT_LUCY 
REAL(kind = double) :: TKIN,PRESSURE,PV_KIN,PV_POT,TAU_TEMP,TAU_TEMP_SQD 
REAL(kind = double) :: PXX,PYY,PXY,PYZ,PZZ,DEGFREE 
REAL(kind = double) :: PKXX,PKYY,PKXY,PKYZ,PKZZ 

INTEGER :: NTS, N_POWER, M_POWER, ISTART, ISUB, MOVWRITE, N_POWER_SUB1, M_POWER_SUB1
INTEGER :: NEQ,N,NPAIRS,MAXPAIR,IPOT_TYPE,SS_POWER
INTEGER :: NCELLX,NCELLY,NCELLZ,NCELL,MAPSIZE,ITHERM 
INTEGER :: ITIME 
 
CHARACTER(7) :: XYZ_FILE_STAT

NAMELIST /input_deck/N, NTS, ISTART,IPOT_TYPE,N_POWER,M_POWER,RCUT, &
                    & TREQ,DENSITY,SHEAR,SIGMA, EPS,SS_POWER,ITHERM,&
                    & DT,ISTART,ISUB,MOVWRITE,TAU_TEMP

END MODULE header_file
