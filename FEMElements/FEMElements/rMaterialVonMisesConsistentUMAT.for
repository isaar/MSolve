C
C Test SUBROUTINE for von Mises material with linear isotropic hardening
C [With Consistent Tangent Constitutive Matrix]
C
C Written by Kostas Papanikolopoulos, 04/11/2004
C
      SUBROUTINE VMCUMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C The above lines are needed by ABAQUS
C User variables are declared below
C
      DIMENSION STR_EL(6),D_EL(6,6),STR_DEV(6)
C STR_EL(6) = Elastic stress predictor
C D_EL(6,6) = Elastic constitutive matrix
C STR_DEV(6)= (New) Deviatoric stress tensor --> Used to form DDSDDE
	DIMENSION ARRAY_K(6,6)
	DATA ARRAY_K   / 0., -1., -1., 0., 0., 0.,
     &				-1.,  0., -1., 0., 0., 0.,
     &				-1., -1.,  0., 0., 0., 0.,
     &				 0.,  0.,  0., .5, 0., 0.,
     &				 0.,  0.,  0., 0., .5, 0.,
     &				 0.,  0.,  0., 0., 0., .5  /
C ARRAY_K(6,6) = Array needed for consistent constitutive matrix formulation
C
C Check the size of the stress or strain component array  
      IF(NTENS.NE.6) THEN
	  WRITE(6,*) "FATAL ERROR IN UMAT. SEE .msg FILE FOR DETAILS"
	  WRITE(7,*) "FATAL ERROR IN UMAT: NTENS IS NOT 6"
	  WRITE(7,*) "PLEASE USE ONLY CONTINUUM (SOLID) ELEMENTS"
	  WRITE(7,*) "WITH THIS MATERIAL MODEL."
	  WRITE(7,*) "EXITING"
c	  CALL XIT
	  STOP
	ENDIF
C
C Evaluate elastic constitutive matrix
C
      DO K1=1,NTENS
        DO K2=1,NTENS
          D_EL(K2,K1) = 0.
        END DO
      END DO
C PROPS(1)=E (Young's modulus), PROPS(2)=nu (Poisson's ratio)
C Do not use nu=0.5 (incrompressible solid)
	C1 = PROPS(1)*(1.-PROPS(2))/(1.+PROPS(2))/(1.-2.*PROPS(2))
	C2 = PROPS(1)*PROPS(2)/(1.+PROPS(2))/(1.-2.*PROPS(2))
	C3 = PROPS(1)/(1.+PROPS(2))/2.
C
      D_EL(1,1) = C1
      D_EL(2,2) = C1
      D_EL(3,3) = C1
C
      D_EL(1,2) = C2
      D_EL(1,3) = C2
      D_EL(2,1) = C2
      D_EL(2,3) = C2
      D_EL(3,1) = C2
      D_EL(3,2) = C2
C
      D_EL(4,4) = C3
      D_EL(5,5) = C3
      D_EL(6,6) = C3
C
C EVALUATE NEW STRESS TENSOR and STATE VARIABLE
C
C Calculate trial (elastic) stress state
	STR_EL = STRESS+MATMUL(D_EL,DSTRAN)
C Yield function <>0 ?
C PROPS(3)=sigma(y) (initial yield stress),
C PROPS(4)=h (hardening modulus for linear hardening)
C STATEV(1)=epsilon(p)_BAR (accumulated plastic strain)
	CJ2_EL = ( (STR_EL(1)-STR_EL(2))**2.+(STR_EL(2)-STR_EL(3))**2.+
     &           (STR_EL(3)-STR_EL(1))**2. )/6.+
     &         STR_EL(4)**2.+STR_EL(5)**2.+STR_EL(6)**2.
	IF (SQRT(3.*CJ2_EL)-(PROPS(3)+PROPS(4)*STATEV(1)).LE.0.) THEN
	  STRESS=STR_EL     ! No plastic correction needed
	  DDSDDE=D_EL       ! Const. matrix = Elastic const. matrix
	  GOTO 900          ! No change to the accumulated plastic strain
	ENDIF
C Return mapping
	G=PROPS(1)/2./(1.+PROPS(2))    ! G = Shear modulus
	DEPSILONP = (SQRT(3.*CJ2_EL)-(PROPS(3)+PROPS(4)*STATEV(1)))/
     &            (3.*G+PROPS(4))
	STATEV(1)=STATEV(1)+DEPSILONP  ! Update accumulated plastic strain
	CCC = G*DEPSILONP*SQRT(3./CJ2_EL)   ! 2. and 1/2. cancel out
	P_EL = (STR_EL(1)+STR_EL(2)+STR_EL(3))/3.  ! P_EL=mean stress (elastic)
	STRESS(1) = STR_EL(1)-CCC*(STR_EL(1)-P_EL)
	STRESS(2) = STR_EL(2)-CCC*(STR_EL(2)-P_EL)
	STRESS(3) = STR_EL(3)-CCC*(STR_EL(3)-P_EL)
	STRESS(4) = STR_EL(4)-CCC*STR_EL(4)
	STRESS(5) = STR_EL(5)-CCC*STR_EL(5)
	STRESS(6) = STR_EL(6)-CCC*STR_EL(6)
C
C CREATE NEW JACOBIAN
C
	P_NEW = (STRESS(1)+STRESS(2)+STRESS(3))/3.
	STR_DEV(1) = STRESS(1)-P_NEW
	STR_DEV(2) = STRESS(2)-P_NEW
	STR_DEV(3) = STRESS(3)-P_NEW
	STR_DEV(4) = STRESS(4)
	STR_DEV(5) = STRESS(5)
	STR_DEV(6) = STRESS(6)
C
	CJ2_NEW = ( (STRESS(1)-STRESS(2))**2.+(STRESS(2)-STRESS(3))**2.+
     &            (STRESS(3)-STRESS(1))**2. )/6.+
     &          STRESS(4)**2.+STRESS(5)**2.+STRESS(6)**2.
	CCC=3.*G*G/(PROPS(4)+3.*G)/CJ2_NEW
	DDD=2.*SQRT(3.)*G*G*DEPSILONP/SQRT(CJ2_EL)
	EEE=.5/CJ2_NEW

      DO K1=1,NTENS
        DO K2=1,NTENS
          DDSDDE(K2,K1) = D_EL(K2,K1)-CCC*STR_DEV(K2)*STR_DEV(K1)
     &				-DDD*(ARRAY_K(K2,K1)-EEE*STR_DEV(K2)*STR_DEV(K1))
        END DO
      END DO
C
  900 CONTINUE
	RETURN
	END