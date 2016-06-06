Subroutine MaterialVonMisesConsistent(afdStrains, fEpsilonP, afStresses, afProps, afStressesNew, fEpsilonPNew, afE)
  !DEC$ ATTRIBUTES DLLEXPORT::MaterialVonMisesConsistent
	Implicit None
    Real, Intent(IN) :: fEpsilonP, afdStrains(6), afStresses(6), afProps(5)
	Real, Intent(OUT) :: afStressesNew(6), fEpsilonPNew, afE(6, 6)
	Real :: fDummy, afDummy(6)
	Character(len=80) :: cName

! Main routine
	afStressesNew = afStresses
	fEpsilonPNew = fEpsilonP
!      SUBROUTINE VMCUMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
!    1 RPL,DDSDDT,DRPLDE,DRPLDT,
!     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
!     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
!     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	Call VMCUMAT(afStressesNew, fEpsilonPNew, afE, fDummy, fDummy, fDummy, &
		fDummy, afDummy, afDummy, fDummy, &
		afDummy, afdStrains, afDummy, fDummy, fDummy, fDummy, afDummy, afDummy, cName, &
		0, 0, 6, 1, afProps, 5, afDummy, fDummy, fDummy, &
		fDummy, afDummy, afDummy, 0, 0, fDummy, fDummy, fDummy, fDummy)

!      CHARACTER*80 CMNAME
!      DIMENSION STRESS(NTENS),STATEV(NSTATV),
!     1 DDSDDE(NTENS,NTENS),
!     2 DDSDDT(NTENS),DRPLDE(NTENS),
!     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
!     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
End Subroutine