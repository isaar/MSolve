Subroutine CalcH20u8pProperties(fPerm, fPorosity, fPoreA, fKs, fKf, fap, faXw, faSw, faQInv, faPermeability)
	Implicit None
	Real, Intent(IN) :: fPerm, fPorosity, fPoreA, fKs, fKf, fap(8)
	Real, Intent(OUT) :: faXw(8), faSw(8), faQInv(8), faPermeability(8)

! Main routine 
	faXw = 1.
	faSw = 1.
	faQInv = 0. + fPorosity * faSw / fKf + (fPoreA - fPorosity) * faXw / fKs
	faPermeability = fPerm
End Subroutine

! ******************************************
! Integer function CalcOnlyKHexa8()
!
! IN:
! iElemNo(integer): Position in the elements matrix 
!					of the element that the stiffness
!					matrix will be calculated.
!
! OUT:
! ---
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate the element stiffness matrix
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20u8pGaussMatrices(iInt, faXYZ, faWeight20, faS20, faDS20, faJ20, faDetJ20, faB20, &
	faDetJ8, faB8, faWeight8, faS8, faDS8)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pGaussMatrices
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faXYZ(3, 20)
	Real, Intent(OUT) :: faDS20(60, iInt ** 3), faS20(20, iInt ** 3), faS8(8, iInt ** 3), &
		faB20(6, 60, iInt ** 3), faB8(6, 24, iInt ** 3), faDetJ20(iInt ** 3), faJ20(3, 3, iInt ** 3), &
		faWeight20(iInt ** 3), faWeight8(iInt ** 3), faDetJ8(iInt ** 3), faDS8(24, iInt ** 3)
	Real fXi, fEta, fZeta, faXYZ8(3, 8), faJ8(3, 3, iInt ** 3), fWt
	Integer M, iX, iY, iZ, iErr, CalcH20JDetJ, CalcH8JDetJ, I, J
	Real GP(4, 4), GW(4, 4)
! GP Matrix stores Gauss - Legendre sampling points
    Data GP/ &
		0.D0,   0.D0,   0.D0,   0.D0,   -.5773502691896D0, &
		.5773502691896D0,   0.D0,   0.D0,   -.7745966692415D0,   0.D0, &
		.7745966692415D0,   0.D0,   -.8611363115941D0, &
		-.3399810435849D0,   .3399810435849D0,   .8611363115941D0 /
! GW Matrix stores Gauss - Legendre weighting factors
	Data GW/ &   
		2.D0,   0.D0,   0.D0,   0.D0,   1.D0,   1.D0, &
		0.D0,  0.D0,    .5555555555556D0,   .8888888888889D0, &
		.5555555555556D0,   0.D0,   .3478548451375D0,   .6521451548625D0, &
		.6521451548625D0, 	 .3478548451375D0 /

! Main routine
	faXYZ8(1, 1) = faXYZ(1, 1)
	faXYZ8(1, 2) = faXYZ(1, 2)
	faXYZ8(1, 3) = faXYZ(1, 3)
	faXYZ8(1, 4) = faXYZ(1, 4)
	faXYZ8(1, 5) = faXYZ(1, 5)
	faXYZ8(1, 6) = faXYZ(1, 6)
	faXYZ8(1, 7) = faXYZ(1, 7)
	faXYZ8(1, 8) = faXYZ(1, 8)
	faXYZ8(2, 1) = faXYZ(2, 1)
	faXYZ8(2, 2) = faXYZ(2, 2)
	faXYZ8(2, 3) = faXYZ(2, 3)
	faXYZ8(2, 4) = faXYZ(2, 4)
	faXYZ8(2, 5) = faXYZ(2, 5)
	faXYZ8(2, 6) = faXYZ(2, 6)
	faXYZ8(2, 7) = faXYZ(2, 7)
	faXYZ8(2, 8) = faXYZ(2, 8)
	faXYZ8(3, 1) = faXYZ(3, 1)
	faXYZ8(3, 2) = faXYZ(3, 2)
	faXYZ8(3, 3) = faXYZ(3, 3)
	faXYZ8(3, 4) = faXYZ(3, 4)
	faXYZ8(3, 5) = faXYZ(3, 5)
	faXYZ8(3, 6) = faXYZ(3, 6)
	faXYZ8(3, 7) = faXYZ(3, 7)
	faXYZ8(3, 8) = faXYZ(3, 8)

	M = 0
	Do iX = 1, iInt
		fXi = GP(iX, iInt)
		Do iY = 1, iInt
			fEta = GP(iY, iInt)
			Do iZ = 1, iInt
				fZeta = GP(iZ, iInt)
				fWt = GW(iX, iInt) * GW(iY, iInt) * GW(iZ, iInt)
				M = M + 1
				Call CalcH20NablaShape(fXi, fEta, fZeta, faDS20(1, M))
				Call CalcH20Shape(fXi, fEta, fZeta, faS20(1, M))
				iErr = CalcH20JDetJ(faXYZ, faDS20(1, M), faJ20(1, 1, M), faDetJ20(M))
				Call CalcH20B(faDetJ20(M), faJ20(1, 1, M), faDS20(1, M), faB20(1, 1, M))
				faWeight20(M) = fWt * faDetJ20(M)

				Call CalcH8NablaShape(fXi, fEta, fZeta, faDS8(1, M))
				Call CalcH8Shape(fXi, fEta, fZeta, faS8(1, M))
				iErr = CalcH8JDetJ(faXYZ8, faDS8(1, M), faJ8(1, 1, M), faDetJ8(M))
				Call CalcH8B(faDetJ8(M), faJ8(1, 1, M), faDS8(1, M), faB8(1, 1, M))
				faWeight8(M) = fWt * faDetJ8(M)
			EndDo
		EndDo
	EndDo
End Subroutine

! ******************************************
! Integer function CalcOnlyMHexa8(iElemNo)
!
! IN:
! iElemNo(integer): Position in the elements matrix 
!					of the element that the stiffness
!					matrix will be calculated.
!
! OUT:
! ---
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate the element stiffness matrix
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20u8pH(iInt, faPermeability, faS, faB, faWeight, faH)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pH
	Implicit None
	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faWeight(iInt ** 3), faS(8, iInt ** 3), faB(6, 24, iInt ** 3), faPermeability(8)
	Real, Intent(OUT) :: faH(36)
	Real :: fStiff, fTerm, fPermeability
	Integer :: I, J, K, M, KS, iX, iY, iZ, iR1, iR2, iR3, iC

! Main routine 
! Calculate Stress-Strain equations
! Calculate stiffness matrix
	M = 0
	Do I = 1, 36
		faH(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
! Add contribution to element stiffness
				fPermeability = 0.
				Do J = 1, 8
					fPermeability = fPermeability + faPermeability(J) * faS(J, M)
				End Do
				KS = 0
				fTerm = faWeight(M) * fPermeability
				Do I = 1, 8
					iR1 = 3 * (I - 1) + 1
					iR2 = iR1 + 1
					iR3 = iR1 + 2
					Do J = I, 8
						KS = KS + 1
						iC = 3 * (J - 1) + 1
						fStiff = faB(1, iC, M) * faB(1, iR1, M) + faB(2, iC + 1, M) * faB(2, iR2, M) + &
							faB(3, iC + 2, M) * faB(3, iR3, M)
						faH(KS) = faH(KS) + fStiff * fTerm
					EndDo
				EndDo
			EndDo
		EndDo
	EndDo
End Subroutine

! ******************************************
! Integer function CalcOnlyMHexa8(iElemNo)
!
! IN:
! iElemNo(integer): Position in the elements matrix 
!					of the element that the stiffness
!					matrix will be calculated.
!
! OUT:
! ---
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate the element stiffness matrix
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20u8pS(iInt, faXwDivQ, faS, faWeight, faM)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pS
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faXwDivQ(8), faWeight(iInt ** 3), faS(8, iInt ** 3)
	Real, Intent(OUT) :: faM(36)
	Integer :: I, J, M, KS, iX, iY, iZ
	Real :: fTerm, fXwDivQ

! Main routine
	M = 0
	Do I = 1, 36
		faM(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
! Add contribution to element stiffness
				fXwDivQ = 0.
				Do J = 1, 8
					fXwDivQ = fXwDivQ + faXwDivQ(J) * faS(J, M)
				End Do
				KS = 0
				fTerm = faWeight(M) * fXwDivQ
				Do I = 1, 8
					Do J = I, 8
						KS = KS + 1
						faM(KS) = faM(KS) + fTerm * faS(I, M) * faS(J, M)
					EndDo
				EndDo
			EndDo
		EndDo
	EndDo

!	Open(100, FILE = 'D:\GAnal\GAnalF\Mass3.txt', STATUS = 'NEW')
!	Do I = 1, 300
!		Write(100, 999) (I), faM(I)
!	EndDo
!	Close(100)
999 FORMAT (I5, F10.5)
End Subroutine

! ******************************************
! Integer function CalcOnlyMHexa8(iElemNo)
!
! IN:
! iElemNo(integer): Position in the elements matrix 
!					of the element that the stiffness
!					matrix will be calculated.
!
! OUT:
! ---
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate the element stiffness matrix
!
! NOTES:
! ---
! ******************************************
!Subroutine CalcH20u8pQTildeQ(iInt, fPoreA, fXw, faB20, faS8, faWeight20, faQ, faQTilde)
Subroutine CalcH20u8pQMinus(iInt, fPoreA, faXw, faB20, faS8, faWeight20, faQTilde)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pQMinus
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fPoreA, faXw(8), faWeight20(iInt ** 3), &
		faB20(6, 60, iInt ** 3), faS8(8, iInt ** 3)
	Real, Intent(OUT) :: faQTilde(60, 8)
	Integer :: I, J, K, M, KS, iX, iY, iZ, iWidth
	Real :: fWeight, fWeightShape, fXw

! Main routine
	M = 0
	Do J = 1, 8
		Do I = 1, 60
			faQTilde(I, J) = 0.
		EndDo
	EndDo

! Compute QTilde
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				fXw = 0.
				Do J = 1, 8
					fXw = fXw + faXw(J) * faS8(J, M)
				End Do
				fWeight = faWeight20(M) * fPoreA * fXw
				Do J = 1, 8
					fWeightShape = fWeight * faS8(J, M)
					Do I = 1, 60
!						faQTilde(I, J) = faQTilde(I, J) + faB20(mod(I - 1, 3) + 1, I, M) * fWeightShape
						faQTilde(I, J) = faQTilde(I, J) - fWeightShape * &
							(faB20(1, I, M) + faB20(2, I, M) + faB20(3, I, M))
					End Do
				End Do
			EndDo
		EndDo
	EndDo

!	faQTilde = fXw * faQTilde
! Compute Q
!	Do J = 1, 8
!		Do I = 1, 60
!			faQ(I, J) = faQTilde(I, J) * fXw
!		EndDo
!	EndDo
!	Open(100, FILE = 'D:\GAnal\GAnalF\Mass3.txt', STATUS = 'NEW')
!	Do I = 1, 300
!		Write(100, 999) (I), faM(I)
!	EndDo
!	Close(100)
999 FORMAT (I5, F10.5)
End Subroutine

! ******************************************
! Integer function CalcOnlyMHexa8(iElemNo)
!
! IN:
! iElemNo(integer): Position in the elements matrix 
!					of the element that the stiffness
!					matrix will be calculated.
!
! OUT:
! ---
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate the element stiffness matrix
!
! NOTES:
! ---
! ******************************************
!Subroutine CalcH20u8pQTildeQ(iInt, fPoreA, fXw, faB20, faS8, faWeight20, faQ, faQTilde)
Subroutine CalcH8u8pQMinus(iInt, fPoreA, faXw, faB20, faS8, faWeight20, faQTilde)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8u8pQMinus
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fPoreA, faXw(8), faWeight20(iInt ** 3), &
		faB20(6, 24, iInt ** 3), faS8(8, iInt ** 3)
	Real, Intent(OUT) :: faQTilde(24, 8)
	Integer :: I, J, K, M, KS, iX, iY, iZ, iWidth
	Real :: fWeight, fWeightShape, fXw

! Main routine
	M = 0
	Do J = 1, 8
		Do I = 1, 24
			faQTilde(I, J) = 0.
		EndDo
	EndDo

! Compute QTilde
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				fXw = 0.
				Do J = 1, 8
					fXw = fXw + faXw(J) * faS8(J, M)
				End Do
				fWeight = faWeight20(M) * fPoreA * fXw
				Do J = 1, 8
					fWeightShape = fWeight * faS8(J, M)
					Do I = 1, 24
!						faQTilde(I, J) = faQTilde(I, J) + faB20(mod(I - 1, 3) + 1, I, M) * fWeightShape
						faQTilde(I, J) = faQTilde(I, J) - fWeightShape * &
							(faB20(1, I, M) + faB20(2, I, M) + faB20(3, I, M))
					End Do
				End Do
			EndDo
		EndDo
	EndDo
999 FORMAT (I5, F10.5)
End Subroutine

Subroutine CalcH20u8pQTilde(iInt, fPoreA, faB20, faS8, faWeight20, faQTilde)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pQTilde
	Implicit None
	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fPoreA, faWeight20(iInt ** 3), &
		faB20(6, 60, iInt ** 3), faS8(8, iInt ** 3)
	Real, Intent(OUT) :: faQTilde(60, 8)
	Integer :: I, J, K, M, KS, iX, iY, iZ, iWidth
	Real :: fWeight, fWeightShape

! Main routine
	M = 0
	Do J = 1, 8
		Do I = 1, 60
			faQTilde(I, J) = 0.
		EndDo
	EndDo

! Compute QTilde
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				fWeight = faWeight20(M) * fPoreA
				Do J = 1, 8
					fWeightShape = fWeight * faS8(J, M)
					Do I = 1, 60
						faQTilde(I, J) = faQTilde(I, J) + fWeightShape * &
							(faB20(1, I, M) + faB20(2, I, M) + faB20(3, I, M))
					End Do
				End Do
			EndDo
		EndDo
	EndDo
End Subroutine

Subroutine CalcH8u8pQTilde(iInt, fPoreA, faB20, faS8, faWeight20, faQTilde)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8u8pQTilde
	Implicit None
	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fPoreA, faWeight20(iInt ** 3), &
		faB20(6, 24, iInt ** 3), faS8(8, iInt ** 3)
	Real, Intent(OUT) :: faQTilde(24, 8)
	Integer :: I, J, K, M, KS, iX, iY, iZ, iWidth
	Real :: fWeight, fWeightShape

! Main routine
	M = 0
	Do J = 1, 8
		Do I = 1, 24
			faQTilde(I, J) = 0.
		EndDo
	EndDo

! Compute QTilde
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				fWeight = faWeight20(M) * fPoreA
				Do J = 1, 8
					fWeightShape = fWeight * faS8(J, M)
					Do I = 1, 24
						faQTilde(I, J) = faQTilde(I, J) + fWeightShape * &
							(faB20(1, I, M) + faB20(2, I, M) + faB20(3, I, M))
					End Do
				End Do
			EndDo
		EndDo
	EndDo
End Subroutine

Subroutine CalcH20u8pForcesFromPoresQ(piEqLM, pfElementQ, pfP, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementQ(60, 8), pfP(8)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iTemp
	Real :: fLocalForces(60)

! Main function
! Initialize local forces vector
	fLocalForces(1) = 0D0
	fLocalForces(2) = 0D0
	fLocalForces(3) = 0D0
	fLocalForces(4) = 0D0
	fLocalForces(5) = 0D0
	fLocalForces(6) = 0D0
	fLocalForces(7) = 0D0
	fLocalForces(8) = 0D0
	fLocalForces(9) = 0D0
	fLocalForces(10) = 0D0
	fLocalForces(11) = 0D0
	fLocalForces(12) = 0D0
	fLocalForces(13) = 0D0
	fLocalForces(14) = 0D0
	fLocalForces(15) = 0D0
	fLocalForces(16) = 0D0
	fLocalForces(17) = 0D0
	fLocalForces(18) = 0D0
	fLocalForces(19) = 0D0
	fLocalForces(20) = 0D0
	fLocalForces(21) = 0D0
	fLocalForces(22) = 0D0
	fLocalForces(23) = 0D0
	fLocalForces(24) = 0D0
	fLocalForces(25) = 0D0
	fLocalForces(26) = 0D0
	fLocalForces(27) = 0D0
	fLocalForces(28) = 0D0
	fLocalForces(29) = 0D0
	fLocalForces(30) = 0D0
	fLocalForces(31) = 0D0
	fLocalForces(32) = 0D0
	fLocalForces(33) = 0D0
	fLocalForces(34) = 0D0
	fLocalForces(35) = 0D0
	fLocalForces(36) = 0D0
	fLocalForces(37) = 0D0
	fLocalForces(38) = 0D0
	fLocalForces(39) = 0D0
	fLocalForces(40) = 0D0
	fLocalForces(41) = 0D0
	fLocalForces(42) = 0D0
	fLocalForces(43) = 0D0
	fLocalForces(44) = 0D0
	fLocalForces(45) = 0D0
	fLocalForces(46) = 0D0
	fLocalForces(47) = 0D0
	fLocalForces(48) = 0D0
	fLocalForces(49) = 0D0
	fLocalForces(50) = 0D0
	fLocalForces(51) = 0D0
	fLocalForces(52) = 0D0
	fLocalForces(53) = 0D0
	fLocalForces(54) = 0D0
	fLocalForces(55) = 0D0
	fLocalForces(56) = 0D0
	fLocalForces(57) = 0D0
	fLocalForces(58) = 0D0
	fLocalForces(59) = 0D0
	fLocalForces(60) = 0D0

! Multiply element Q with vector pfP and store to fLocalForces
	Do J = 1, 8
		Do I = 1, 60
			fLocalForces(I) = fLocalForces(I) + pfElementQ(I, J) * pfP(J)
		End Do
	End Do

! Update the global loads vector (pfGlobalForces)
	iTemp = piEqLM(1)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(1)
	iTemp = piEqLM(2)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(2)
	iTemp = piEqLM(3)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(3)
	iTemp = piEqLM(4)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(4)
	iTemp = piEqLM(5)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(5)
	iTemp = piEqLM(6)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(6)
	iTemp = piEqLM(7)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(7)
	iTemp = piEqLM(8)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(8)
	iTemp = piEqLM(9)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(9)
	iTemp = piEqLM(10)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(10)
	iTemp = piEqLM(11)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(11)
	iTemp = piEqLM(12)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(12)
	iTemp = piEqLM(13)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(13)
	iTemp = piEqLM(14)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(14)
	iTemp = piEqLM(15)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(15)
	iTemp = piEqLM(16)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(16)
	iTemp = piEqLM(17)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(17)
	iTemp = piEqLM(18)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(18)
	iTemp = piEqLM(19)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(19)
	iTemp = piEqLM(20)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(20)
	iTemp = piEqLM(21)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(21)
	iTemp = piEqLM(22)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(22)
	iTemp = piEqLM(23)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(23)
	iTemp = piEqLM(24)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(24)
	iTemp = piEqLM(25)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(25)
	iTemp = piEqLM(26)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(26)
	iTemp = piEqLM(27)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(27)
	iTemp = piEqLM(28)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(28)
	iTemp = piEqLM(29)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(29)
	iTemp = piEqLM(30)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(30)
	iTemp = piEqLM(31)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(31)
	iTemp = piEqLM(32)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(32)
	iTemp = piEqLM(33)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(33)
	iTemp = piEqLM(34)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(34)
	iTemp = piEqLM(35)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(35)
	iTemp = piEqLM(36)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(36)
	iTemp = piEqLM(37)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(37)
	iTemp = piEqLM(38)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(38)
	iTemp = piEqLM(39)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(39)
	iTemp = piEqLM(40)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(40)
	iTemp = piEqLM(41)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(41)
	iTemp = piEqLM(42)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(42)
	iTemp = piEqLM(43)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(43)
	iTemp = piEqLM(44)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(44)
	iTemp = piEqLM(45)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(45)
	iTemp = piEqLM(46)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(46)
	iTemp = piEqLM(47)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(47)
	iTemp = piEqLM(48)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(48)
	iTemp = piEqLM(49)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(49)
	iTemp = piEqLM(50)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(50)
	iTemp = piEqLM(51)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(51)
	iTemp = piEqLM(52)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(52)
	iTemp = piEqLM(53)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(53)
	iTemp = piEqLM(54)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(54)
	iTemp = piEqLM(55)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(55)
	iTemp = piEqLM(56)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(56)
	iTemp = piEqLM(57)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(57)
	iTemp = piEqLM(58)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(58)
	iTemp = piEqLM(59)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(59)
	iTemp = piEqLM(60)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(60)
End Subroutine

Subroutine CalcH20u8pForcesFromPoresH(piEqLM, pfElementH, pfP, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementH(36), pfP(8)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iRow, iAccR, iTemp, iDataPos
	Real :: fLocalForces(8)

! Main function
! Initialize local forces vector
	fLocalForces(1) = 0D0
	fLocalForces(2) = 0D0
	fLocalForces(3) = 0D0
	fLocalForces(4) = 0D0
	fLocalForces(5) = 0D0
	fLocalForces(6) = 0D0
	fLocalForces(7) = 0D0
	fLocalForces(8) = 0D0

! Multiply SKYLINE Element H with vector pfP and store to fLocalForces
	iDataPos = 1
	Do I = 1, 8
		fLocalForces(I) = fLocalForces(I) + pfElementH(iDataPos) * pfP(I)
		iDataPos = iDataPos + 1
		Do J = I + 1, 8
			fLocalForces(I) = fLocalForces(I) + pfElementH(iDataPos) * pfP(J)
			fLocalForces(J) = fLocalForces(J) + pfElementH(iDataPos) * pfP(I)
			iDataPos = iDataPos + 1
		End Do
	End Do

! Update the global loads vector (pfGlobalForces)
	iTemp = piEqLM(1)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(1)
	iTemp = piEqLM(2)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(2)
	iTemp = piEqLM(3)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(3)
	iTemp = piEqLM(4)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(4)
	iTemp = piEqLM(5)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(5)
	iTemp = piEqLM(6)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(6)
	iTemp = piEqLM(7)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(7)
	iTemp = piEqLM(8)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(8)
End Subroutine

Subroutine CalcH20u8pForcesFromAccWater(iInt, iImpermeable, piImpermeableIDs, piEqLM, &
	piEqIDs, fDPS, pfElementB, pfAcc, pfWeights, pfGlobalForces)
	Implicit None
	Integer, Intent(IN) :: iInt, iImpermeable, piEqLM(1), piEqIDs(4, 1), piImpermeableIDs(1)
	Real, Intent(IN) :: fDPS, pfElementB(1), pfAcc(3), pfWeights(1)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, iX, iY, iZ, iRow, iAccR, iTemp, iDataPos, iWPos, piEqs(8)
	Real :: fLocalForces(8), fAccDPS(3), fWeight, fTerm

! Main function
! Initialize local forces vector
	fLocalForces(1) = 0D0
	fLocalForces(2) = 0D0
	fLocalForces(3) = 0D0
	fLocalForces(4) = 0D0
	fLocalForces(5) = 0D0
	fLocalForces(6) = 0D0
	fLocalForces(7) = 0D0
	fLocalForces(8) = 0D0
	piEqs(1) = piEqLM(1)
	piEqs(2) = piEqLM(2)
	piEqs(3) = piEqLM(3)
	piEqs(4) = piEqLM(4)
	piEqs(5) = piEqLM(5)
	piEqs(6) = piEqLM(6)
	piEqs(7) = piEqLM(7)
	piEqs(8) = piEqLM(8)

! Multiply body forces with saturation, permeability and density
	fAccDPS(1) = fDPS * pfAcc(1)
	fAccDPS(2) = fDPS * pfAcc(2)
	fAccDPS(3) = fDPS * pfAcc(3)

! Multiply Element B with vector pfAcc and store to fLocalForces
	iDataPos = 1
	iWPos = 0
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				iWPos = iWPos + 1
				fWeight = pfWeights(iWPos)
				Do I = 1, 8
					fTerm = pfElementB(iDataPos) * fAccDPS(1) + &
						pfElementB(iDataPos + 1) * fAccDPS(2) + &
						pfElementB(iDataPos + 2) * fAccDPS(3)
					fLocalForces(I) = fLocalForces(I) + fTerm * fWeight
!					fLocalForces(I) = fLocalForces(I) + pfElementB(iDataPos) * fAccDPS(1) * &
!						fTemp
					iDataPos = iDataPos + 3
				End Do
			End Do
		End Do
	End Do

! Update the global loads vector (pfGlobalForces)
	Do I = 1, iImpermeable
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(1)) piEqs(1) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(2)) piEqs(2) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(3)) piEqs(3) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(4)) piEqs(4) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(5)) piEqs(5) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(6)) piEqs(6) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(7)) piEqs(7) = 0
!		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(8)) piEqs(8) = 0
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(1)) fLocalForces(1) = fLocalForces(1) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(2)) fLocalForces(2) = fLocalForces(2) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(3)) fLocalForces(3) = fLocalForces(3) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(4)) fLocalForces(4) = fLocalForces(4) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(5)) fLocalForces(5) = fLocalForces(5) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(6)) fLocalForces(6) = fLocalForces(6) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(7)) fLocalForces(7) = fLocalForces(7) / 1D11
		if (piEqIDs(4, piImpermeableIDs(I)).EQ.piEqs(8)) fLocalForces(8) = fLocalForces(8) / 1D11
	End Do

	iTemp = piEqs(1)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(1)
	iTemp = piEqs(2)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(2)
	iTemp = piEqs(3)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(3)
	iTemp = piEqs(4)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(4)
	iTemp = piEqs(5)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(5)
	iTemp = piEqs(6)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(6)
	iTemp = piEqs(7)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(7)
	iTemp = piEqs(8)
	If (iTemp.GT.0) pfGlobalForces(iTemp) = pfGlobalForces(iTemp) + fLocalForces(8)
End Subroutine

!DEC$ ATTRIBUTES DLLEXPORT::CalcH20u8pForcesWaterAcc
Subroutine CalcH20u8pForcesWaterAcc(iInt, alImpermeable, ffDensity, afPermeability, afXw, afSw, afAcc, afS, afB, &
	afWeights, afLocalForces)
	Implicit None
	Integer, Intent(IN) :: iInt
	Logical, Intent(IN) :: alImpermeable(8)
	Real, Intent(IN) :: ffDensity, afPermeability(8), afXw(8), afSw(8), afS(8, iInt ** 3), &
		afB(6, 24, iInt ** 3), afAcc(3), afWeights(*)
	Real, Intent(OUT) :: afLocalForces(8)
	Integer :: I, iX, iY, iZ, iDataPos, iWPos
	Real :: afAccDPS(3), fWeight, fTerm

! Main function
! Initialize local forces vector
	afLocalForces(1) = 0D0
	afLocalForces(2) = 0D0
	afLocalForces(3) = 0D0
	afLocalForces(4) = 0D0
	afLocalForces(5) = 0D0
	afLocalForces(6) = 0D0
	afLocalForces(7) = 0D0
	afLocalForces(8) = 0D0

! Multiply Element B with vector pfAcc and store to fLocalForces
	iDataPos = 1
	iWPos = 0
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				iWPos = iWPos + 1
				fTerm = 0.
				Do I = 1, 8
					fTerm = fTerm + afPermeability(I) * afXw(I) * afSw(I) * afS(I, iWPos)
				End Do
				afAccDPS = afAcc * ffDensity * fTerm
				fWeight = afWeights(iWPos)
				Do I = 0, 7
					fTerm = afB(1, (I * 3) + 1, iWPos) * afAccDPS(1) + &
						afB(2, (I * 3) + 2, iWPos) * afAccDPS(2) + &
						afB(3, (I * 3) + 3, iWPos) * afAccDPS(3)
!					fTerm = afElementB(iDataPos) * afAccDPS(1) + &
!						afElementB(iDataPos + 1) * afAccDPS(2) + &
!						afElementB(iDataPos + 2) * afAccDPS(3)
					afLocalForces(I + 1) = afLocalForces(I + 1) + fTerm * fWeight
!					fLocalForces(I) = fLocalForces(I) + pfElementB(iDataPos) * fAccDPS(1) * fTemp
					iDataPos = iDataPos + 3
				End Do
			End Do
		End Do
	End Do

! If impermeable, nullify force
	Do I = 1, 8
		If (alImpermeable(I)) afLocalForces(I) = 0.
	End Do
End Subroutine

