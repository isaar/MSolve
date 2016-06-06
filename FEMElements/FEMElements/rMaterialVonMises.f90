!	afStresses:	Current afStresses (IN), New afStresses (OUT)
!	fEpsilonP:	Current State variables (IN), New State variables (OUT)
!	afE:	Constitutive Matrix (OUT)
!	afdStrains:	INCREMENTAL Strains (IN)
!	afProps:	Properties (IN)
! afStressesElastic(6) = Elastic afStresses predictor
! afDElastic(6,6) = Elastic constitutive matrix
! afStressesDeviatoric(6)= (New) Deviatoric Stresses tensor --> Used to form afE

Subroutine MaterialVonMises(afdStrains, fEpsilonP, afStresses, afProps, afStressesNew, fEpsilonPNew, afE)
  !DEC$ ATTRIBUTES DLLEXPORT::MaterialVonMises
	Implicit None
    Real, Intent(IN) :: fEpsilonP, afdStrains(6), afStresses(6), afProps(5)
	Real, Intent(OUT) :: afStressesNew(6), fEpsilonPNew, afE(6, 6)

	Real :: afStressesElastic(6), afDElastic(6, 6), afStressesDeviatoric(6)
	Integer :: K1, K2
	Real :: fTemp1, fTemp2, fTemp3, fG, fdEpsilonP, P_EL, P_NEW, CJ2_EL, CJ2_NEW

! afProps(1)=E (Young's modulus), afProps(2)=nu (Poisson's ratio)
! Do not use nu=0.5 (incrompressible solid)
! Evaluate elastic constitutive matrix
    afDElastic = 0.
	fTemp1 = afProps(1) * (1. - afProps(2)) / (1. + afProps(2)) / (1. - 2. * afProps(2))
	fTemp2 = afProps(1) * afProps(2) / (1. + afProps(2)) / (1. - 2. * afProps(2))
	fTemp3 = afProps(1) / (1. + afProps(2)) / 2.
	afDElastic(1, 1) = fTemp1
	afDElastic(2, 2) = fTemp1
	afDElastic(3, 3) = fTemp1
	afDElastic(1, 2) = fTemp2
	afDElastic(1, 3) = fTemp2
	afDElastic(2, 1) = fTemp2
	afDElastic(2, 3) = fTemp2
	afDElastic(3, 1) = fTemp2
	afDElastic(3, 2) = fTemp2
	afDElastic(4, 4) = fTemp3
	afDElastic(5, 5) = fTemp3
	afDElastic(6, 6) = fTemp3

! EVALUATE NEW afStresses TENSOR and STATE VARIABLE
! Calculate trial (elastic) afStresses state
	afStressesElastic = afStresses + MatMul(afDElastic, afdStrains)
! Yield function <>0 ?
! afProps(3)=sigma(y) (initial yield afStresses),
! afProps(4)=h (hardening modulus for linear hardening)
! fEpsilonP(1)=epsilon(p)_BAR (accumulated plastic strain)
	CJ2_EL = ((afStressesElastic(1) - afStressesElastic(2)) ** 2 + (afStressesElastic(2) - afStressesElastic(3)) ** 2 + &
		(afStressesElastic(3) - afStressesElastic(1)) ** 2) / 6. + afStressesElastic(4) ** 2 + afStressesElastic(5) ** 2 + &
		afStressesElastic(6) ** 2
	fTemp1 = SQRT(3. * CJ2_EL) - (afProps(3) + afProps(4) * fEpsilonP)
	If (fTemp1 <= 0.) Then
		afStressesNew = afStressesElastic     ! No plastic correction needed
		afE = afDElastic       ! Const. matrix = Elastic const. matrix, no change to the accumulated plastic strain
	Else
! Return mapping
		fG = afProps(1) / 2. / (1. + afProps(2))    ! fG = Shear modulus
		fdEpsilonP = fTemp1 / (3. * fG + afProps(4))
		fEpsilonPNew = fEpsilonP + fdEpsilonP  ! Update accumulated plastic strain
		fTemp1 = fG * fdEpsilonP * SQRT(3. / CJ2_EL)   ! 2. and 1/2. cancel out
		P_EL = (afStressesElastic(1) + afStressesElastic(2) + afStressesElastic(3)) / 3.  ! P_EL=mean afStresses (elastic)
		afStressesNew(1) = afStressesElastic(1) - fTemp1 * (afStressesElastic(1) - P_EL)
		afStressesNew(2) = afStressesElastic(2) - fTemp1 * (afStressesElastic(2) - P_EL)
		afStressesNew(3) = afStressesElastic(3) - fTemp1 * (afStressesElastic(3) - P_EL)
		afStressesNew(4) = afStressesElastic(4) - fTemp1 * afStressesElastic(4)
		afStressesNew(5) = afStressesElastic(5) - fTemp1 * afStressesElastic(5)
		afStressesNew(6) = afStressesElastic(6) - fTemp1 * afStressesElastic(6)

! CREATE NEW JACOBIAN
		P_NEW = (afStressesNew(1) + afStressesNew(2) + afStressesNew(3)) / 3.
		afStressesDeviatoric(1) = afStressesNew(1) - P_NEW
		afStressesDeviatoric(2) = afStressesNew(2) - P_NEW
		afStressesDeviatoric(3) = afStressesNew(3) - P_NEW
		afStressesDeviatoric(4) = afStressesNew(4)
		afStressesDeviatoric(5) = afStressesNew(5)
		afStressesDeviatoric(6) = afStressesNew(6)

		CJ2_NEW = ((afStressesNew(1) - afStressesNew(2)) ** 2 + (afStressesNew(2) - afStressesNew(3)) ** 2 + &
			(afStressesNew(3) - afStressesNew(1)) ** 2 ) / 6. + afStressesNew(4) ** 2 + afStressesNew(5) ** 2 + &
			afStressesNew(6) ** 2
		fTemp1 = 3. * fG * fG / (afProps(4) + 3. * fG) / CJ2_NEW
		Do K1 = 1, 6
			Do K2 = 1, 6
				afE(K2, K1) = afDElastic(K2, K1) - fTemp1 * afStressesDeviatoric(K2) * afStressesDeviatoric(K1)
			End Do
		End Do
!		afE = afDElastic       ! Const. matrix = Elastic const. matrix, no change to the accumulated plastic strain
	End If
End Subroutine
