! ******************************************
! Integer function CalcH8NablaShape()
!
! IN:
! iElemType(integer): Element type (see hCommon.h)
! XYZ(float), XI(float), ETA(float), ZETA(float), 
! B(float), DET(float)
!
! OUT:
! nDOFs(integer): Number of DOFs per element
! nProps(itneger): Number of properties per property type
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
!    Evaluate the derivatives of the shape functions
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20NablaShape(fXi, fEta, fZeta, faDS)
	Implicit None

	Real, Intent(IN) :: fXi, fEta, fZeta
	Real, Intent(OUT) :: faDS(60)
	Real fXiP, fEtaP, fZetaP, fXiM, fEtaM, fZetaM

! Main routine
	fXiP = (1D0 + fXi)
	fEtaP = (1D0 + fEta)
	fZetaP = (1D0 + fZeta)
	fXiM = (1D0 - fXi)
	fEtaM = (1D0 - fEta)
	fZetaM = (1D0 - fZeta)

! Calculation of natural coordinate derivatives of the shape functions
! Corresponding to fXi
	faDS(1) = 0.125D0 * fEtaM * fZetaM * (2D0 * fXi + fEta + fZeta + 1D0) 
	faDS(2) = 0.125D0 * fEtaM * fZetaM * (2D0 * fXi - fEta - fZeta - 1D0) 
	faDS(3) = 0.125D0 * fEtaP * fZetaM * (2D0 * fXi + fEta - fZeta - 1D0) 
	faDS(4) = 0.125D0 * fEtaP * fZetaM * (2D0 * fXi - fEta + fZeta + 1D0) 
	faDS(5) = 0.125D0 * fEtaM * fZetaP * (2D0 * fXi + fEta - fZeta + 1D0) 
	faDS(6) = 0.125D0 * fEtaM * fZetaP * (2D0 * fXi - fEta + fZeta - 1D0) 
	faDS(7) = 0.125D0 * fEtaP * fZetaP * (2D0 * fXi + fEta + fZeta - 1D0) 
	faDS(8) = 0.125D0 * fEtaP * fZetaP * (2D0 * fXi - fEta - fZeta + 1D0) 
	faDS(9) = -0.5D0 * fXi * fEtaM * fZetaM
	faDS(10) = 0.25D0 * fEtaM * fEtaP * fZetaM 
	faDS(11) = -0.5D0 * fXi * fEtaP * fZetaM 
	faDS(12) = -0.25D0 * fEtaM * fEtaP * fZetaM 
	faDS(13) = -0.5D0 * fXi * fEtaM * fZetaP 
	faDS(14) = 0.25D0 * fEtaM * fEtaP * fZetaP 
	faDS(15) = -0.5D0 * fXi * fEtaP * fZetaP 
	faDS(16) = -0.25D0 * fEtaM * fEtaP * fZetaP 
	faDS(17) = -0.25D0 * fEtaM * fZetaP * fZetaM 
	faDS(18) = 0.25D0 * fEtaM * fZetaP * fZetaM
	faDS(19) = 0.25D0 * fEtaP * fZetaP * fZetaM 
	faDS(20) = -0.25D0 * fEtaP * fZetaP * fZetaM

! Corresponding to fEta
	faDS(21) = 0.125D0 * fXiM * fZetaM * (fXi + 2D0 * fEta + fZeta + 1D0) 
	faDS(22) = 0.125D0 * fXiP * fZetaM * (-fXi + 2D0 * fEta + fZeta + 1D0) 
	faDS(23) = 0.125D0 * fXiP * fZetaM * (fXi + 2D0 * fEta - fZeta - 1D0) 
	faDS(24) = 0.125D0 * fXiM * fZetaM * (-fXi + 2D0 * fEta - fZeta - 1D0) 
	faDS(25) = 0.125D0 * fXiM * fZetaP * (fXi + 2D0 * fEta - fZeta + 1D0) 
	faDS(26) = 0.125D0 * fXiP * fZetaP * (-fXi + 2D0 * fEta - fZeta + 1D0) 
	faDS(27) = 0.125D0 * fXiP * fZetaP * (fXi + 2D0 * fEta + fZeta - 1D0) 
	faDS(28) = 0.125D0 * fXiM * fZetaP * (-fXi + 2D0 * fEta + fZeta - 1D0) 
	faDS(29) = -0.25D0 * fXiP * fXiM * fZetaM
	faDS(30) = -0.5D0 * fXiP * fEta * fZetaM 
	faDS(31) = 0.25D0 * fXiM * fXiP * fZetaM 
	faDS(32) = -0.5D0 * fXiM * fEta * fZetaM 
	faDS(33) = -0.25D0 * fXiP * fXiM * fZetaP 
	faDS(34) = -0.5D0 * fEta * fXiP * fZetaP 
	faDS(35) = 0.25D0 * fXiM * fXiP * fZetaP 
	faDS(36) = -0.5D0 * fXiM * fEta * fZetaP 
	faDS(37) = -0.25D0 * fXiM * fZetaP * fZetaM 
	faDS(38) = -0.25D0 * fXiP * fZetaP * fZetaM 
	faDS(39) = 0.25D0 * fXiP * fZetaP * fZetaM 
	faDS(40) = 0.25D0 * fXiM * fZetaP * fZetaM

! Corresponding to fZeta
	faDS(41) = 0.125D0 * fXiM * fEtaM * (fXi + fEta + 2D0 * fZeta + 1D0) 
	faDS(42) = 0.125D0 * fXiP * fEtaM * (-fXi + fEta + 2D0 * fZeta + 1D0) 
	faDS(43) = 0.125D0 * fXiP * fEtaP * (-fXi - fEta + 2D0 * fZeta + 1D0) 
	faDS(44) = 0.125D0 * fXiM * fEtaP * (fXi - fEta + 2D0 * fZeta + 1D0) 
	faDS(45) = 0.125D0 * fXiM * fEtaM * (-fXi - fEta + 2D0 * fZeta - 1D0) 
	faDS(46) = 0.125D0 * fXiP * fEtaM * (fXi - fEta + 2D0 * fZeta - 1D0) 
	faDS(47) = 0.125D0 * fXiP * fEtaP * (fXi + fEta + 2D0 * fZeta - 1D0) 
	faDS(48) = 0.125D0 * fXiM * fEtaP * (-fXi + fEta + 2D0 * fZeta - 1D0) 
	faDS(49) = -0.25D0 * fXiP * fXiM * fEtaM 
	faDS(50) = -0.25D0 * fXiP * fEtaP * fEtaM 
	faDS(51) = -0.25D0 * fXiM * fXiP * fEtaP 
	faDS(52) = -0.25D0 * fXiM * fEtaP * fEtaM
	faDS(53) = 0.25D0 * fXiP * fXiM * fEtaM 
	faDS(54) = 0.25D0 * fXiP * fEtaP * fEtaM 
	faDS(55) = 0.25D0 * fXiM * fXiP * fEtaP 
	faDS(56) = 0.25D0 * fXiM * fEtaM * fEtaP 
	faDS(57) = -0.5D0 * fXiM * fEtaM * fZeta 
	faDS(58) = -0.5D0 * fXiP * fEtaM * fZeta 
	faDS(59) = -0.5D0 * fXiP * fEtaP * fZeta 
	faDS(60) = -0.5D0 * fXiM * fEtaP * fZeta
End Subroutine

! ******************************************
! Integer function CalcShapeHexa8()
!
! IN:
! iElemType(integer): Element type (see hCommon.h)
!
! OUT:
! nDOFs(integer): Number of DOFs per element
! nProps(itneger): Number of properties per property type
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate shape functions for calculation of consistent mass matrix
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20Shape(fXi, fEta, fZeta, faS)
	Implicit None

	Real, Intent(IN) :: fXi, fEta, fZeta
	Real, Intent(OUT) :: faS(20)
	Real, Parameter :: fSqC125 = 0.5D0
	Real fXiP, fEtaP, fZetaP, fXiM, fEtaM, fZetaM, &
		fXi2, fEta2, fZeta2

! Main routine
	fXiP = (1D0 + fXi) * fSqC125
	fEtaP = (1D0 + fEta) * fSqC125
	fZetaP = (1D0 + fZeta) * fSqC125
	fXiM = (1D0 - fXi) * fSqC125
	fEtaM = (1D0 - fEta) * fSqC125
	fZetaM = (1D0 - fZeta) * fSqC125
	fXi2 = (1D0 - fXi * fXi)
	fEta2 = (1D0 - fEta * fEta)
	fZeta2 = (1D0 - fZeta * fZeta)

	faS(1) = fXiM * fEtaM * fZetaM * (-fXi - fEta - fZeta - 2D0)
	faS(2) = fXiP * fEtaM * fZetaM * (fXi - fEta - fZeta - 2D0)
	faS(3) = fXiP * fEtaP * fZetaM * (fXi + fEta - fZeta - 2D0)
	faS(4) = fXiM * fEtaP * fZetaM * (-fXi + fEta - fZeta - 2D0)
	faS(5) = fXiM * fEtaM * fZetaP * (-fXi - fEta + fZeta - 2D0)
	faS(6) = fXiP * fEtaM * fZetaP * (fXi - fEta + fZeta - 2D0)
	faS(7) = fXiP * fEtaP * fZetaP * (fXi + fEta + fZeta - 2D0)
	faS(8) = fXiM * fEtaP * fZetaP * (-fXi + fEta + fZeta - 2D0)
	faS(9) = fXi2 * fEtaM * fZetaM
	faS(10) = fEta2 * fXiP * fZetaM
	faS(11) = fXi2 * fEtaP * fZetaM
	faS(12) = fEta2 * fXiM * fZetaM
	faS(13) = fXi2 * fEtaM * fZetaP
	faS(14) = fEta2 * fXiP * fZetaP
	faS(15) = fXi2 * fEtaP * fZetaP
	faS(16) = fEta2 * fXiM * fZetaP
	faS(17) = fZeta2 * fXiM * fEtaM
	faS(18) = fZeta2 * fXiP * fEtaM
	faS(19) = fZeta2 * fXiP * fEtaP
	faS(20) = fZeta2 * fXiM * fEtaP
End Subroutine

! ******************************************
! Integer function CalcShapeHexa8()
!
! IN:
! iElemType(integer): Element type (see hCommon.h)
!
! OUT:
! nDOFs(integer): Number of DOFs per element
! nProps(itneger): Number of properties per property type
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate shape functions for calculation of consistent mass matrix
!
! NOTES:
! ---
! ******************************************
Integer Function CalcH20JDetJ(faXYZ, faDS, faJ, fDetJ)
	Implicit None
	Include 'hCommon.h'

	Real, Intent(IN) :: faXYZ(3, 20), faDS(60)
	Real, Intent(OUT) :: faJ(3, 3), fDetJ
	Real fDet1, fDet2, fDet3

! Main routine
	CalcH20JDetJ = iErrNegDet
! Compute the Jacobian at (XI, ETA, ZETA)
	faJ(1, 1) = faDS(1) * faXYZ(1, 1) + &
		faDS(2) * faXYZ(1, 2) + faDS(3) * faXYZ(1, 3) + &
		faDS(4) * faXYZ(1, 4) + faDS(5) * faXYZ(1, 5) + &
		faDS(6) * faXYZ(1, 6) + faDS(7) * faXYZ(1, 7) + &
		faDS(8) * faXYZ(1, 8) + faDS(9) * faXYZ(1, 9) + &
		faDS(10) * faXYZ(1, 10) + faDS(11) * faXYZ(1, 11) + &
		faDS(12) * faXYZ(1, 12) + faDS(13) * faXYZ(1, 13) + &
		faDS(14) * faXYZ(1, 14) + faDS(15) * faXYZ(1, 15) + &
		faDS(16) * faXYZ(1, 16) + faDS(17) * faXYZ(1, 17) + &
		faDS(18) * faXYZ(1, 18) + faDS(19) * faXYZ(1, 19) + &
		faDS(20) * faXYZ(1, 20)
	faJ(1, 2) = faDS(1) * faXYZ(2, 1) + &
		faDS(2) * faXYZ(2, 2) + faDS(3) * faXYZ(2, 3) + &
		faDS(4) * faXYZ(2, 4) + faDS(5) * faXYZ(2, 5) + &
		faDS(6) * faXYZ(2, 6) + faDS(7) * faXYZ(2, 7) + &
		faDS(8) * faXYZ(2, 8) + faDS(9) * faXYZ(2, 9) + &
		faDS(10) * faXYZ(2, 10) + faDS(11) * faXYZ(2, 11) + &
		faDS(12) * faXYZ(2, 12) + faDS(13) * faXYZ(2, 13) + &
		faDS(14) * faXYZ(2, 14) + faDS(15) * faXYZ(2, 15) + &
		faDS(16) * faXYZ(2, 16) + faDS(17) * faXYZ(2, 17) + &
		faDS(18) * faXYZ(2, 18) + faDS(19) * faXYZ(2, 19) + &
		faDS(20) * faXYZ(2, 20)
	faJ(1, 3) = faDS(1) * faXYZ(3, 1) + &
		faDS(2) * faXYZ(3, 2) + faDS(3) * faXYZ(3, 3) + &
		faDS(4) * faXYZ(3, 4) + faDS(5) * faXYZ(3, 5) + &
		faDS(6) * faXYZ(3, 6) + faDS(7) * faXYZ(3, 7) + &
		faDS(8) * faXYZ(3, 8) + faDS(9) * faXYZ(3, 9) + &
		faDS(10) * faXYZ(3, 10) + faDS(11) * faXYZ(3, 11) + &
		faDS(12) * faXYZ(3, 12) + faDS(13) * faXYZ(3, 13) + &
		faDS(14) * faXYZ(3, 14) + faDS(15) * faXYZ(3, 15) + &
		faDS(16) * faXYZ(3, 16) + faDS(17) * faXYZ(3, 17) + &
		faDS(18) * faXYZ(3, 18) + faDS(19) * faXYZ(3, 19) + &
		faDS(20) * faXYZ(3, 20)
	faJ(2, 1) = faDS(21) * faXYZ(1, 1) + &
		faDS(22) * faXYZ(1, 2) + faDS(23) * faXYZ(1, 3) + &
		faDS(24) * faXYZ(1, 4) + faDS(25) * faXYZ(1, 5) + &
		faDS(26) * faXYZ(1, 6) + faDS(27) * faXYZ(1, 7) + &
		faDS(28) * faXYZ(1, 8) + faDS(29) * faXYZ(1, 9) + &
		faDS(30) * faXYZ(1, 10) + faDS(31) * faXYZ(1, 11) + &
		faDS(32) * faXYZ(1, 12) + faDS(33) * faXYZ(1, 13) + &
		faDS(34) * faXYZ(1, 14) + faDS(35) * faXYZ(1, 15) + &
		faDS(36) * faXYZ(1, 16) + faDS(37) * faXYZ(1, 17) + &
		faDS(38) * faXYZ(1, 18) + faDS(39) * faXYZ(1, 19) + &
		faDS(40) * faXYZ(1, 20)
	faJ(2, 2) = faDS(21) * faXYZ(2, 1) + &
		faDS(22) * faXYZ(2, 2) + faDS(23) * faXYZ(2, 3) + &
		faDS(24) * faXYZ(2, 4) + faDS(25) * faXYZ(2, 5) + &
		faDS(26) * faXYZ(2, 6) + faDS(27) * faXYZ(2, 7) + &
		faDS(28) * faXYZ(2, 8) + faDS(29) * faXYZ(2, 9) + &
		faDS(30) * faXYZ(2, 10) + faDS(31) * faXYZ(2, 11) + &
		faDS(32) * faXYZ(2, 12) + faDS(33) * faXYZ(2, 13) + &
		faDS(34) * faXYZ(2, 14) + faDS(35) * faXYZ(2, 15) + &
		faDS(36) * faXYZ(2, 16) + faDS(37) * faXYZ(2, 17) + &
		faDS(38) * faXYZ(2, 18) + faDS(39) * faXYZ(2, 19) + &
		faDS(40) * faXYZ(2, 20)
	faJ(2, 3) = faDS(21) * faXYZ(3, 1) + &
		faDS(22) * faXYZ(3, 2) + faDS(23) * faXYZ(3, 3) + &
		faDS(24) * faXYZ(3, 4) + faDS(25) * faXYZ(3, 5) + &
		faDS(26) * faXYZ(3, 6) + faDS(27) * faXYZ(3, 7) + &
		faDS(28) * faXYZ(3, 8) + faDS(29) * faXYZ(3, 9) + &
		faDS(30) * faXYZ(3, 10) + faDS(31) * faXYZ(3, 11) + &
		faDS(32) * faXYZ(3, 12) + faDS(33) * faXYZ(3, 13) + &
		faDS(34) * faXYZ(3, 14) + faDS(35) * faXYZ(3, 15) + &
		faDS(36) * faXYZ(3, 16) + faDS(37) * faXYZ(3, 17) + &
		faDS(38) * faXYZ(3, 18) + faDS(39) * faXYZ(3, 19) + &
		faDS(40) * faXYZ(3, 20)
	faJ(3, 1) = faDS(41) * faXYZ(1, 1) + &
		faDS(42) * faXYZ(1, 2) + faDS(43) * faXYZ(1, 3) + &
		faDS(44) * faXYZ(1, 4) + faDS(45) * faXYZ(1, 5) + &
		faDS(46) * faXYZ(1, 6) + faDS(47) * faXYZ(1, 7) + &
		faDS(48) * faXYZ(1, 8) + faDS(49) * faXYZ(1, 9) + &
		faDS(50) * faXYZ(1, 10) + faDS(51) * faXYZ(1, 11) + &
		faDS(52) * faXYZ(1, 12) + faDS(53) * faXYZ(1, 13) + &
		faDS(54) * faXYZ(1, 14) + faDS(55) * faXYZ(1, 15) + &
		faDS(56) * faXYZ(1, 16) + faDS(57) * faXYZ(1, 17) + &
		faDS(58) * faXYZ(1, 18) + faDS(59) * faXYZ(1, 19) + &
		faDS(60) * faXYZ(1, 20)
	faJ(3, 2) = faDS(41) * faXYZ(2, 1) + &
		faDS(42) * faXYZ(2, 2) + faDS(43) * faXYZ(2, 3) + &
		faDS(44) * faXYZ(2, 4) + faDS(45) * faXYZ(2, 5) + &
		faDS(46) * faXYZ(2, 6) + faDS(47) * faXYZ(2, 7) + &
		faDS(48) * faXYZ(2, 8) + faDS(49) * faXYZ(2, 9) + &
		faDS(50) * faXYZ(2, 10) + faDS(51) * faXYZ(2, 11) + &
		faDS(52) * faXYZ(2, 12) + faDS(53) * faXYZ(2, 13) + &
		faDS(54) * faXYZ(2, 14) + faDS(55) * faXYZ(2, 15) + &
		faDS(56) * faXYZ(2, 16) + faDS(57) * faXYZ(2, 17) + &
		faDS(58) * faXYZ(2, 18) + faDS(59) * faXYZ(2, 19) + &
		faDS(60) * faXYZ(2, 20)
	faJ(3, 3) = faDS(41) * faXYZ(3, 1) + &
		faDS(42) * faXYZ(3, 2) + faDS(43) * faXYZ(3, 3) + &
		faDS(44) * faXYZ(3, 4) + faDS(45) * faXYZ(3, 5) + &
		faDS(46) * faXYZ(3, 6) + faDS(47) * faXYZ(3, 7) + &
		faDS(48) * faXYZ(3, 8) + faDS(49) * faXYZ(3, 9) + &
		faDS(50) * faXYZ(3, 10) + faDS(51) * faXYZ(3, 11) + &
		faDS(52) * faXYZ(3, 12) + faDS(53) * faXYZ(3, 13) + &
		faDS(54) * faXYZ(3, 14) + faDS(55) * faXYZ(3, 15) + &
		faDS(56) * faXYZ(3, 16) + faDS(57) * faXYZ(3, 17) + &
		faDS(58) * faXYZ(3, 18) + faDS(59) * faXYZ(3, 19) + &
		faDS(60) * faXYZ(3, 20)
	fDet1 = faJ(1, 1) * (faJ(2, 2) * faJ(3, 3) - faJ(3, 2) * faJ(2, 3))
	fDet2 = -faJ(1, 2) * (faJ(2, 1) * faJ(3, 3) - faJ(3, 1) * faJ(2, 3))
	fDet3 = faJ(1, 3) * (faJ(2, 1) * faJ(3, 2) - faJ(3, 1) * faJ(2, 2))
	fDetJ = fDet1 + fDet2 + fDet3
! Error if determinant is negative (return)
	If (fDetJ.LT.0.00000001) then
		Write(*, *) "Negative determinant (H20)"
		Return
	EndIf
	CalcH20JDetJ = iErrNoError
End Function

! ******************************************
! Integer function CalcShapeHexa8()
!
! IN:
! iElemType(integer): Element type (see hCommon.h)
!
! OUT:
! nDOFs(integer): Number of DOFs per element
! nProps(itneger): Number of properties per property type
!
! RETURNS:
! (integer): Number of nodes per element
!
! UPDATES:
! ---
!
! PURPOSE:
! Calculate shape functions for calculation of consistent mass matrix
!
! NOTES:
! ---
! ******************************************
Subroutine CalcH20B(fDetJ, faJ, faDS, faB)
	Implicit None
	Include 'hCommon.h'

	Real, Intent(IN) :: fDetJ, faJ(3, 3), faDS(60)
	Real, Intent(OUT) :: faB(6, 60)
	Real faJInv(3, 3), fDetInv

! Main routine
! Calculate the Jacobian Inverse (faJInv)
	fDetInv = 1D0 / fDetJ
	faJInv(1, 1) = (faJ(2, 2) * faJ(3, 3) - faJ(3, 2) * faJ(2, 3)) * fDetInv
	faJInv(2, 1) = (faJ(3, 1) * faJ(2, 3) - faJ(2, 1) * faJ(3, 3)) * fDetInv
	faJInv(3, 1) = (faJ(2, 1) * faJ(3, 2) - faJ(3, 1) * faJ(2, 2)) * fDetInv
	faJInv(1, 2) = (faJ(3, 2) * faJ(1, 3) - faJ(1, 2) * faJ(3, 3)) * fDetInv
	faJInv(2, 2) = (faJ(1, 1) * faJ(3, 3) - faJ(3, 1) * faJ(1, 3)) * fDetInv
	faJInv(3, 2) = (faJ(3, 1) * faJ(1, 2) - faJ(3, 2) * faJ(1, 1)) * fDetInv
	faJInv(1, 3) = (faJ(1, 2) * faJ(2, 3) - faJ(2, 2) * faJ(1, 3)) * fDetInv
	faJInv(2, 3) = (faJ(2, 1) * faJ(1, 3) - faJ(1, 1) * faJ(2, 3)) * fDetInv
	faJInv(3, 3) = (faJ(1, 1) * faJ(2, 2) - faJ(2, 1) * faJ(1, 2)) * fDetInv

! Calculate the global derivative operator [faB]
	faB(1, 1) = faJInv(1, 1) * faDS(1) + faJInv(1, 2) * faDS(21) + faJInv(1, 3) * faDS(41)
	faB(1, 2) = 0D0
	faB(1, 3) = 0D0
	faB(1, 4) = faJInv(1, 1) * faDS(2) + faJInv(1, 2) * faDS(22) + faJInv(1, 3) * faDS(42)
	faB(1, 5) = 0D0
	faB(1, 6) = 0D0
	faB(1, 7) = faJInv(1, 1) * faDS(3) + faJInv(1, 2) * faDS(23) + faJInv(1, 3) * faDS(43)
	faB(1, 8) = 0D0
	faB(1, 9) = 0D0
	faB(1, 10) = faJInv(1, 1) * faDS(4) + faJInv(1, 2) * faDS(24) + faJInv(1, 3) * faDS(44)
	faB(1, 11) = 0D0
	faB(1, 12) = 0D0
	faB(1, 13) = faJInv(1, 1) * faDS(5) + faJInv(1, 2) * faDS(25) + faJInv(1, 3) * faDS(45)
	faB(1, 14) = 0D0
	faB(1, 15) = 0D0
	faB(1, 16) = faJInv(1, 1) * faDS(6) + faJInv(1, 2) * faDS(26) + faJInv(1, 3) * faDS(46)
	faB(1, 17) = 0D0
	faB(1, 18) = 0D0
	faB(1, 19) = faJInv(1, 1) * faDS(7) + faJInv(1, 2) * faDS(27) + faJInv(1, 3) * faDS(47)
	faB(1, 20) = 0D0
	faB(1, 21) = 0D0
	faB(1, 22) = faJInv(1, 1) * faDS(8) + faJInv(1, 2) * faDS(28) + faJInv(1, 3) * faDS(48)
	faB(1, 23) = 0D0
	faB(1, 24) = 0D0
	faB(1, 25) = faJInv(1, 1) * faDS(9) + faJInv(1, 2) * faDS(29) + faJInv(1, 3) * faDS(49)
	faB(1, 26) = 0D0
	faB(1, 27) = 0D0
	faB(1, 28) = faJInv(1, 1) * faDS(10) + faJInv(1, 2) * faDS(30) + faJInv(1, 3) * faDS(50)
	faB(1, 29) = 0D0
	faB(1, 30) = 0D0
	faB(1, 31) = faJInv(1, 1) * faDS(11) + faJInv(1, 2) * faDS(31) + faJInv(1, 3) * faDS(51)
	faB(1, 32) = 0D0
	faB(1, 33) = 0D0
	faB(1, 34) = faJInv(1, 1) * faDS(12) + faJInv(1, 2) * faDS(32) + faJInv(1, 3) * faDS(52)
	faB(1, 35) = 0D0
	faB(1, 36) = 0D0
	faB(1, 37) = faJInv(1, 1) * faDS(13) + faJInv(1, 2) * faDS(33) + faJInv(1, 3) * faDS(53)
	faB(1, 38) = 0D0
	faB(1, 39) = 0D0
	faB(1, 40) = faJInv(1, 1) * faDS(14) + faJInv(1, 2) * faDS(34) + faJInv(1, 3) * faDS(54)
	faB(1, 41) = 0D0
	faB(1, 42) = 0D0
	faB(1, 43) = faJInv(1, 1) * faDS(15) + faJInv(1, 2) * faDS(35) + faJInv(1, 3) * faDS(55)
	faB(1, 44) = 0D0
	faB(1, 45) = 0D0
	faB(1, 46) = faJInv(1, 1) * faDS(16) + faJInv(1, 2) * faDS(36) + faJInv(1, 3) * faDS(56)
	faB(1, 47) = 0D0
	faB(1, 48) = 0D0
	faB(1, 49) = faJInv(1, 1) * faDS(17) + faJInv(1, 2) * faDS(37) + faJInv(1, 3) * faDS(57)
	faB(1, 50) = 0D0
	faB(1, 51) = 0D0
	faB(1, 52) = faJInv(1, 1) * faDS(18) + faJInv(1, 2) * faDS(38) + faJInv(1, 3) * faDS(58)
	faB(1, 53) = 0D0
	faB(1, 54) = 0D0
	faB(1, 55) = faJInv(1, 1) * faDS(19) + faJInv(1, 2) * faDS(39) + faJInv(1, 3) * faDS(59)
	faB(1, 56) = 0D0
	faB(1, 57) = 0D0
	faB(1, 58) = faJInv(1, 1) * faDS(20) + faJInv(1, 2) * faDS(40) + faJInv(1, 3) * faDS(60)
	faB(1, 59) = 0D0
	faB(1, 60) = 0D0
	faB(2, 1) = 0D0
	faB(2, 2) = faJInv(2, 1) * faDS(1) + faJInv(2, 2) * faDS(21) + faJInv(2, 3) * faDS(41)
	faB(2, 3) = 0D0
	faB(2, 4) = 0D0
	faB(2, 5) = faJInv(2, 1) * faDS(2) + faJInv(2, 2) * faDS(22) + faJInv(2, 3) * faDS(42)
	faB(2, 6) = 0D0
	faB(2, 7) = 0D0
	faB(2, 8) = faJInv(2, 1) * faDS(3) + faJInv(2, 2) * faDS(23) + faJInv(2, 3) * faDS(43)
	faB(2, 9) = 0D0
	faB(2, 10) = 0D0
	faB(2, 11) = faJInv(2, 1) * faDS(4) + faJInv(2, 2) * faDS(24) + faJInv(2, 3) * faDS(44)
	faB(2, 12) = 0D0
	faB(2, 13) = 0D0
	faB(2, 14) = faJInv(2, 1) * faDS(5) + faJInv(2, 2) * faDS(25) + faJInv(2, 3) * faDS(45)
	faB(2, 15) = 0D0
	faB(2, 16) = 0D0
	faB(2, 17) = faJInv(2, 1) * faDS(6) + faJInv(2, 2) * faDS(26) + faJInv(2, 3) * faDS(46)
	faB(2, 18) = 0D0
	faB(2, 19) = 0D0
	faB(2, 20) = faJInv(2, 1) * faDS(7) + faJInv(2, 2) * faDS(27) + faJInv(2, 3) * faDS(47)
	faB(2, 21) = 0D0
	faB(2, 22) = 0D0
	faB(2, 23) = faJInv(2, 1) * faDS(8) + faJInv(2, 2) * faDS(28) + faJInv(2, 3) * faDS(48)
	faB(2, 24) = 0D0
	faB(2, 25) = 0D0
	faB(2, 26) = faJInv(2, 1) * faDS(9) + faJInv(2, 2) * faDS(29) + faJInv(2, 3) * faDS(49)
	faB(2, 27) = 0D0
	faB(2, 28) = 0D0
	faB(2, 29) = faJInv(2, 1) * faDS(10) + faJInv(2, 2) * faDS(30) + faJInv(2, 3) * faDS(50)
	faB(2, 30) = 0D0
	faB(2, 31) = 0D0
	faB(2, 32) = faJInv(2, 1) * faDS(11) + faJInv(2, 2) * faDS(31) + faJInv(2, 3) * faDS(51)
	faB(2, 33) = 0D0
	faB(2, 34) = 0D0
	faB(2, 35) = faJInv(2, 1) * faDS(12) + faJInv(2, 2) * faDS(32) + faJInv(2, 3) * faDS(52)
	faB(2, 36) = 0D0
	faB(2, 37) = 0D0
	faB(2, 38) = faJInv(2, 1) * faDS(13) + faJInv(2, 2) * faDS(33) + faJInv(2, 3) * faDS(53)
	faB(2, 39) = 0D0
	faB(2, 40) = 0D0
	faB(2, 41) = faJInv(2, 1) * faDS(14) + faJInv(2, 2) * faDS(34) + faJInv(2, 3) * faDS(54)
	faB(2, 42) = 0D0
	faB(2, 43) = 0D0
	faB(2, 44) = faJInv(2, 1) * faDS(15) + faJInv(2, 2) * faDS(35) + faJInv(2, 3) * faDS(55)
	faB(2, 45) = 0D0
	faB(2, 46) = 0D0
	faB(2, 47) = faJInv(2, 1) * faDS(16) + faJInv(2, 2) * faDS(36) + faJInv(2, 3) * faDS(56)
	faB(2, 48) = 0D0
	faB(2, 49) = 0D0
	faB(2, 50) = faJInv(2, 1) * faDS(17) + faJInv(2, 2) * faDS(37) + faJInv(2, 3) * faDS(57)
	faB(2, 51) = 0D0
	faB(2, 52) = 0D0
	faB(2, 53) = faJInv(2, 1) * faDS(18) + faJInv(2, 2) * faDS(38) + faJInv(2, 3) * faDS(58)
	faB(2, 54) = 0D0
	faB(2, 55) = 0D0
	faB(2, 56) = faJInv(2, 1) * faDS(19) + faJInv(2, 2) * faDS(39) + faJInv(2, 3) * faDS(59)
	faB(2, 57) = 0D0
	faB(2, 58) = 0D0
	faB(2, 59) = faJInv(2, 1) * faDS(20) + faJInv(2, 2) * faDS(40) + faJInv(2, 3) * faDS(60)
	faB(2, 60) = 0D0
	faB(3, 1) = 0D0
	faB(3, 2) = 0D0
	faB(3, 3) = faJInv(3, 1) * faDS(1) + faJInv(3, 2) * faDS(21) + faJInv(3, 3) * faDS(41)
	faB(3, 4) = 0D0
	faB(3, 5) = 0D0
	faB(3, 6) = faJInv(3, 1) * faDS(2) + faJInv(3, 2) * faDS(22) + faJInv(3, 3) * faDS(42)
	faB(3, 7) = 0D0
	faB(3, 8) = 0D0
	faB(3, 9) = faJInv(3, 1) * faDS(3) + faJInv(3, 2) * faDS(23) + faJInv(3, 3) * faDS(43)
	faB(3, 10) = 0D0
	faB(3, 11) = 0D0
	faB(3, 12) = faJInv(3, 1) * faDS(4) + faJInv(3, 2) * faDS(24) + faJInv(3, 3) * faDS(44)
	faB(3, 13) = 0D0
	faB(3, 14) = 0D0
	faB(3, 15) = faJInv(3, 1) * faDS(5) + faJInv(3, 2) * faDS(25) + faJInv(3, 3) * faDS(45)
	faB(3, 16) = 0D0
	faB(3, 17) = 0D0
	faB(3, 18) = faJInv(3, 1) * faDS(6) + faJInv(3, 2) * faDS(26) + faJInv(3, 3) * faDS(46)
	faB(3, 19) = 0D0
	faB(3, 20) = 0D0
	faB(3, 21) = faJInv(3, 1) * faDS(7) + faJInv(3, 2) * faDS(27) + faJInv(3, 3) * faDS(47)
	faB(3, 22) = 0D0
	faB(3, 23) = 0D0
	faB(3, 24) = faJInv(3, 1) * faDS(8) + faJInv(3, 2) * faDS(28) + faJInv(3, 3) * faDS(48)
	faB(3, 25) = 0D0
	faB(3, 26) = 0D0
	faB(3, 27) = faJInv(3, 1) * faDS(9) + faJInv(3, 2) * faDS(29) + faJInv(3, 3) * faDS(49)
	faB(3, 28) = 0D0
	faB(3, 29) = 0D0
	faB(3, 30) = faJInv(3, 1) * faDS(10) + faJInv(3, 2) * faDS(30) + faJInv(3, 3) * faDS(50)
	faB(3, 31) = 0D0
	faB(3, 32) = 0D0
	faB(3, 33) = faJInv(3, 1) * faDS(11) + faJInv(3, 2) * faDS(31) + faJInv(3, 3) * faDS(51)
	faB(3, 34) = 0D0
	faB(3, 35) = 0D0
	faB(3, 36) = faJInv(3, 1) * faDS(12) + faJInv(3, 2) * faDS(32) + faJInv(3, 3) * faDS(52)
	faB(3, 37) = 0D0
	faB(3, 38) = 0D0
	faB(3, 39) = faJInv(3, 1) * faDS(13) + faJInv(3, 2) * faDS(33) + faJInv(3, 3) * faDS(53)
	faB(3, 40) = 0D0
	faB(3, 41) = 0D0
	faB(3, 42) = faJInv(3, 1) * faDS(14) + faJInv(3, 2) * faDS(34) + faJInv(3, 3) * faDS(54)
	faB(3, 43) = 0D0
	faB(3, 44) = 0D0
	faB(3, 45) = faJInv(3, 1) * faDS(15) + faJInv(3, 2) * faDS(35) + faJInv(3, 3) * faDS(55)
	faB(3, 46) = 0D0
	faB(3, 47) = 0D0
	faB(3, 48) = faJInv(3, 1) * faDS(16) + faJInv(3, 2) * faDS(36) + faJInv(3, 3) * faDS(56)
	faB(3, 49) = 0D0
	faB(3, 50) = 0D0
	faB(3, 51) = faJInv(3, 1) * faDS(17) + faJInv(3, 2) * faDS(37) + faJInv(3, 3) * faDS(57)
	faB(3, 52) = 0D0
	faB(3, 53) = 0D0
	faB(3, 54) = faJInv(3, 1) * faDS(18) + faJInv(3, 2) * faDS(38) + faJInv(3, 3) * faDS(58)
	faB(3, 55) = 0D0
	faB(3, 56) = 0D0
	faB(3, 57) = faJInv(3, 1) * faDS(19) + faJInv(3, 2) * faDS(39) + faJInv(3, 3) * faDS(59)
	faB(3, 58) = 0D0
	faB(3, 59) = 0D0
	faB(3, 60) = faJInv(3, 1) * faDS(20) + faJInv(3, 2) * faDS(40) + faJInv(3, 3) * faDS(60)
	faB(4, 1) = faB(2, 2)
	faB(4, 2) = faB(1, 1)
	faB(4, 3) = 0D0
	faB(4, 4) = faB(2, 5)
	faB(4, 5) = faB(1, 4)
	faB(4, 6) = 0D0
	faB(4, 7) = faB(2, 8)
	faB(4, 8) = faB(1, 7)
	faB(4, 9) = 0D0
	faB(4, 10) = faB(2, 11)
	faB(4, 11) = faB(1, 10)
	faB(4, 12) = 0D0
	faB(4, 13) = faB(2, 14)
	faB(4, 14) = faB(1, 13)
	faB(4, 15) = 0D0
	faB(4, 16) = faB(2, 17)
	faB(4, 17) = faB(1, 16)
	faB(4, 18) = 0D0
	faB(4, 19) = faB(2, 20)
	faB(4, 20) = faB(1, 19)
	faB(4, 21) = 0D0
	faB(4, 22) = faB(2, 23)
	faB(4, 23) = faB(1, 22)
	faB(4, 24) = 0D0
	faB(4, 25) = faB(2, 26)
	faB(4, 26) = faB(1, 25)
	faB(4, 27) = 0D0
	faB(4, 28) = faB(2, 29)
	faB(4, 29) = faB(1, 28)
	faB(4, 30) = 0D0
	faB(4, 31) = faB(2, 32)
	faB(4, 32) = faB(1, 31)
	faB(4, 33) = 0D0
	faB(4, 34) = faB(2, 35)
	faB(4, 35) = faB(1, 34)
	faB(4, 36) = 0D0
	faB(4, 37) = faB(2, 38)
	faB(4, 38) = faB(1, 37)
	faB(4, 39) = 0D0
	faB(4, 40) = faB(2, 41)
	faB(4, 41) = faB(1, 40)
	faB(4, 42) = 0D0
	faB(4, 43) = faB(2, 44)
	faB(4, 44) = faB(1, 43)
	faB(4, 45) = 0D0
	faB(4, 46) = faB(2, 47)
	faB(4, 47) = faB(1, 46)
	faB(4, 48) = 0D0
	faB(4, 49) = faB(2, 50)
	faB(4, 50) = faB(1, 49)
	faB(4, 51) = 0D0
	faB(4, 52) = faB(2, 53)
	faB(4, 53) = faB(1, 52)
	faB(4, 54) = 0D0
	faB(4, 55) = faB(2, 56)
	faB(4, 56) = faB(1, 55)
	faB(4, 57) = 0D0
	faB(4, 58) = faB(2, 59)
	faB(4, 59) = faB(1, 58)
	faB(4, 60) = 0D0
	faB(5, 1) = 0D0
	faB(5, 2) = faB(3, 3)
	faB(5, 3) = faB(2, 2)
	faB(5, 4) = 0D0
	faB(5, 5) = faB(3, 6)
	faB(5, 6) = faB(2, 5)
	faB(5, 7) = 0D0
	faB(5, 8) = faB(3, 9)
	faB(5, 9) = faB(2, 8)
	faB(5, 10) = 0D0
	faB(5, 11) = faB(3, 12)
	faB(5, 12) = faB(2, 11)
	faB(5, 13) = 0D0
	faB(5, 14) = faB(3, 15)
	faB(5, 15) = faB(2, 14)
	faB(5, 16) = 0D0
	faB(5, 17) = faB(3, 18)
	faB(5, 18) = faB(2, 17)
	faB(5, 19) = 0D0
	faB(5, 20) = faB(3, 21)
	faB(5, 21) = faB(2, 20)
	faB(5, 22) = 0D0
	faB(5, 23) = faB(3, 24)
	faB(5, 24) = faB(2, 23)
	faB(5, 25) = 0D0
	faB(5, 26) = faB(3, 27)
	faB(5, 27) = faB(2, 26)
	faB(5, 28) = 0D0
	faB(5, 29) = faB(3, 30)
	faB(5, 30) = faB(2, 29)
	faB(5, 31) = 0D0
	faB(5, 32) = faB(3, 33)
	faB(5, 33) = faB(2, 32)
	faB(5, 34) = 0D0
	faB(5, 35) = faB(3, 36)
	faB(5, 36) = faB(2, 35)
	faB(5, 37) = 0D0
	faB(5, 38) = faB(3, 39)
	faB(5, 39) = faB(2, 38)
	faB(5, 40) = 0D0
	faB(5, 41) = faB(3, 42)
	faB(5, 42) = faB(2, 41)
	faB(5, 43) = 0D0
	faB(5, 44) = faB(3, 45)
	faB(5, 45) = faB(2, 44)
	faB(5, 46) = 0D0
	faB(5, 47) = faB(3, 48)
	faB(5, 48) = faB(2, 47)
	faB(5, 49) = 0D0
	faB(5, 50) = faB(3, 51)
	faB(5, 51) = faB(2, 50)
	faB(5, 52) = 0D0
	faB(5, 53) = faB(3, 54)
	faB(5, 54) = faB(2, 53)
	faB(5, 55) = 0D0
	faB(5, 56) = faB(3, 57)
	faB(5, 57) = faB(2, 56)
	faB(5, 58) = 0D0
	faB(5, 59) = faB(3, 60)
	faB(5, 60) = faB(2, 59)
	faB(6, 1) = faB(3, 3)
	faB(6, 2) = 0D0
	faB(6, 3) = faB(1, 1)
	faB(6, 4) = faB(3, 6)
	faB(6, 5) = 0D0
	faB(6, 6) = faB(1, 4)
	faB(6, 7) = faB(3, 9)
	faB(6, 8) = 0D0
	faB(6, 9) = faB(1, 7)
	faB(6, 10) = faB(3, 12)
	faB(6, 11) = 0D0
	faB(6, 12) = faB(1, 10)
	faB(6, 13) = faB(3, 15)
	faB(6, 14) = 0D0
	faB(6, 15) = faB(1, 13)
	faB(6, 16) = faB(3, 18)
	faB(6, 17) = 0D0
	faB(6, 18) = faB(1, 16)
	faB(6, 19) = faB(3, 21)
	faB(6, 20) = 0D0
	faB(6, 21) = faB(1, 19)
	faB(6, 22) = faB(3, 24)
	faB(6, 23) = 0D0
	faB(6, 24) = faB(1, 22)
	faB(6, 25) = faB(3, 27)
	faB(6, 26) = 0D0
	faB(6, 27) = faB(1, 25)
	faB(6, 28) = faB(3, 30)
	faB(6, 29) = 0D0
	faB(6, 30) = faB(1, 28)
	faB(6, 31) = faB(3, 33)
	faB(6, 32) = 0D0
	faB(6, 33) = faB(1, 31)
	faB(6, 34) = faB(3, 36)
	faB(6, 35) = 0D0
	faB(6, 36) = faB(1, 34)
	faB(6, 37) = faB(3, 39)
	faB(6, 38) = 0D0
	faB(6, 39) = faB(1, 37)
	faB(6, 40) = faB(3, 42)
	faB(6, 41) = 0D0
	faB(6, 42) = faB(1, 40)
	faB(6, 43) = faB(3, 45)
	faB(6, 44) = 0D0
	faB(6, 45) = faB(1, 43)
	faB(6, 46) = faB(3, 48)
	faB(6, 47) = 0D0
	faB(6, 48) = faB(1, 46)
	faB(6, 49) = faB(3, 51)
	faB(6, 50) = 0D0
	faB(6, 51) = faB(1, 49)
	faB(6, 52) = faB(3, 54)
	faB(6, 53) = 0D0
	faB(6, 54) = faB(1, 52)
	faB(6, 55) = faB(3, 57)
	faB(6, 56) = 0D0
	faB(6, 57) = faB(1, 55)
	faB(6, 58) = faB(3, 60)
	faB(6, 59) = 0D0
	faB(6, 60) = faB(1, 58)
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
Subroutine CalcH20GaussMatrices(iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20GaussMatrices
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faXYZ(3, 20)
	Real, Intent(OUT) :: faDS(60, iInt ** 3), faS(20, iInt ** 3), &
		faB(6, 60, iInt ** 3), faDetJ(iInt ** 3), faJ(3, 3, iInt ** 3), &
		faWeight(iInt ** 3)
	Real fXi, fEta, fZeta
	Integer M, iX, iY, iZ, iErr, CalcH20JDetJ
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
	M = 0
	Do iX = 1, iInt
		fXi = GP(iX, iInt)
		Do iY = 1, iInt
			fEta = GP(iY, iInt)
			Do iZ = 1, iInt
				fZeta = GP(iZ, iInt)
				M = M + 1
				Call CalcH20NablaShape(fXi, fEta, fZeta, faDS(1, M))
				Call CalcH20Shape(fXi, fEta, fZeta, faS(1, M))
				iErr = CalcH20JDetJ(faXYZ, faDS(1, M), faJ(1, 1, M), faDetJ(M))
				Call CalcH20B(faDetJ(M), faJ(1, 1, M), faDS(1, M), faB(1, 1, M))
				faWeight(M) = GW(iX, iInt) * GW(iY, iInt) * GW(iZ, iInt) * faDetJ(M)
			EndDo
		EndDo
	EndDo
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
Subroutine CalcH20K(iInt, faE, faB, faWeight, faK)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20K
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faWeight(iInt ** 3), faB(6, 60, iInt ** 3), faE(6, 6, iInt ** 3)
	Real, Intent(OUT) :: faK(1830)
	Real :: faEB(6), fStiff
	Integer :: I, J, K, M, KS, iX, iY, iZ

! Main routine 
! Calculate stiffness matrix
	M = 0
	Do I = 1, 1830
		faK(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
! Add contribution to element stiffness
				KS = 0
				Do I = 1, 60
					faEB(1) = faE(1, 1, M) * faB(1, I, M) + faE(1, 2, M) * faB(2, I, M) + &
						faE(1, 3, M) * faB(3, I, M) + faE(1, 4, M) * faB(4, I, M) + &
						faE(1, 5, M) * faB(5, I, M) + faE(1, 6, M) * faB(6, I, M)
					faEB(2) = faE(2, 1, M) * faB(1, I, M) + faE(2, 2, M) * faB(2, I, M) + &
						faE(2, 3, M) * faB(3, I, M) + faE(2, 4, M) * faB(4, I, M) + &
						faE(2, 5, M) * faB(5, I, M) + faE(2, 6, M) * faB(6, I, M)
					faEB(3) = faE(3, 1, M) * faB(1, I, M) + faE(3, 2, M) * faB(2, I, M) + &
						faE(3, 3, M) * faB(3, I, M) + faE(3, 4, M) * faB(4, I, M) + &
						faE(3, 5, M) * faB(5, I, M) + faE(3, 6, M) * faB(6, I, M)
					faEB(4) = faE(4, 1, M) * faB(1, I, M) + faE(4, 2, M) * faB(2, I, M) + &
						faE(4, 3, M) * faB(3, I, M) + faE(4, 4, M) * faB(4, I, M) + &
						faE(4, 5, M) * faB(5, I, M) + faE(4, 6, M) * faB(6, I, M)
					faEB(5) = faE(5, 1, M) * faB(1, I, M) + faE(5, 2, M) * faB(2, I, M) + &
						faE(5, 3, M) * faB(3, I, M) + faE(5, 4, M) * faB(4, I, M) + &
						faE(5, 5, M) * faB(5, I, M) + faE(5, 6, M) * faB(6, I, M)
					faEB(6) = faE(6, 1, M) * faB(1, I, M) + faE(6, 2, M) * faB(2, I, M) + &
						faE(6, 3, M) * faB(3, I, M) + faE(6, 4, M) * faB(4, I, M) + &
						faE(6, 5, M) * faB(5, I, M) + faE(6, 6, M) * faB(6, I, M)
					Do J = I, 60
						KS = KS + 1
						fStiff = faB(1, J, M) * faEB(1) + faB(2, J, M) * faEB(2) + &
							faB(3, J, M) * faEB(3) + faB(4, J, M) * faEB(4) + &
							faB(5, J, M) * faEB(5) + faB(6, J, M) * faEB(6)
						faK(KS) = faK(KS) + fStiff * faWeight(M)
					EndDo
				EndDo
			EndDo
		EndDo
	EndDo
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
Subroutine CalcH20Strains(iInt, faB, fau, faStrains)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20Strains
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faB(6, 60, iInt ** 3), fau(60)
	Real, Intent(OUT) :: faStrains(6, iInt ** 3)
	Integer :: M

! Main routine 
	Do M = 1, iInt ** 3
		faStrains(:, M) = MatMul(faB(:, :, M), fau)
	EndDo
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
! AKYRO!!!! POULO H ROUTINA!!!
Subroutine CalcH20Stresses(iInt, faE, faStrains, faStresses)
	Implicit None
	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faE(6, 6), faStrains(6, iInt ** 3)
	Real, Intent(OUT) :: faStresses(6, iInt ** 3)
	Integer :: M

! Main routine 
	Do M = 1, iInt ** 3
		faStresses(:, M) = MatMul(faE, faStrains(:, M))
	EndDo
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
Subroutine CalcH20Forces(iInt, faB, faWeight, faStresses, faForces)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20Forces
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faWeight(iInt ** 3), faB(6, 60, iInt ** 3), faStresses(6, iInt ** 3)
	Real, Intent(OUT) :: faForces(60)
	Integer :: M, iX, iY, iZ

! Main routine 
! Calculate stiffness matrix
	M = 0
	Do iX = 1, 60
		faForces(iX) = 0.
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				faForces = faForces + faWeight(M) * MatMul(Transpose(faB(:, :, M)), faStresses(:, M))
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
Integer Function CalcH20MConsistent(iInt, fDensity, faS, faWeight, faM)
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fDensity, faWeight(iInt ** 3), faS(20, iInt ** 3)
	Real, Intent(OUT) :: faM(1830)
	Integer I, J, K, M, KS, iX, iY, iZ, iWidth
	Real fWeight

! Main routine
	M = 0
	Do I = 1, 1830
		faM(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				KS = -2
				fWeight = faWeight(M) * fDensity
				Do I = 1, 20
					Do J = I, 20
						KS = KS + 3
						faM(KS) = faM(KS) + faS(I, M) * faS(J, M) * fWeight
					EndDo

					Do K = 2, 3
!						iWidth = (21 - I) * 3 + 2 - K
						iWidth = 65 - K - 3 * I
						Do J = I, 20
							KS = KS + 3
							faM(KS) = faM(KS - iWidth)
						EndDo
						KS = KS - K + 1
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
End Function

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
Integer Function CalcH20MLumped(iInt, fDensity, faWeight, faM)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH20MLumped
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fDensity, faWeight(iInt ** 3)
	Real, Intent(OUT) :: faM(1830)
	Integer M, iS2, iS1
	Real fMassTerm

! Main routine
	Do M = 1, 1830
		faM(M) = 0D0
	EndDo

	fMassTerm = 0D0
	Do M = 1, iInt ** 3
		fMassTerm = fMassTerm + faWeight(M)
	EndDo	
	fMassTerm = fMassTerm * fDensity * 0.05D0

	iS2 = 1
	Do iS1 = 60, 1, -1
		faM(iS2) = fMassTerm
		iS2 = iS2 + iS1
	EndDo
End Function

Subroutine CalcH20ForcesFromDisp(piEqLM, pfElementK, pfU, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementK(1), pfU(60)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iRow, iAccR, iTemp, iDataPos
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

! Multiply SKYLINE Element M with vector pfAcc and store to fLocalForces
	iDataPos = 1
	Do I = 1, 60
		fLocalForces(I) = fLocalForces(I) + pfElementK(iDataPos) * pfU(I)
		iDataPos = iDataPos + 1
		Do J = I + 1, 60
			fLocalForces(I) = fLocalForces(I) + pfElementK(iDataPos) * pfU(J)
			fLocalForces(J) = fLocalForces(J) + pfElementK(iDataPos) * pfU(I)
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

Subroutine CalcH20ForcesAcc(iInt, bIsLumped, fDensity, afAcc, afS, afWeights, afLocalForces)
	Implicit None
	Integer, Intent(IN) :: iInt
	Logical, Intent(IN) :: bIsLumped
	Real, Intent(IN) :: fDensity, afS(20, iInt ** 3), afAcc(3), afWeights(*)
	Real, Intent(OUT) :: afLocalForces(60)
	Integer :: I, J, iAccR, iDataPos, CalcH20MLumped, CalcH20MConsistent
	Real :: afM(1830)

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
	afLocalForces(9) = 0D0
	afLocalForces(10) = 0D0
	afLocalForces(11) = 0D0
	afLocalForces(12) = 0D0
	afLocalForces(13) = 0D0
	afLocalForces(14) = 0D0
	afLocalForces(15) = 0D0
	afLocalForces(16) = 0D0
	afLocalForces(17) = 0D0
	afLocalForces(18) = 0D0
	afLocalForces(19) = 0D0
	afLocalForces(20) = 0D0
	afLocalForces(21) = 0D0
	afLocalForces(22) = 0D0
	afLocalForces(23) = 0D0
	afLocalForces(24) = 0D0
	afLocalForces(25) = 0D0
	afLocalForces(26) = 0D0
	afLocalForces(27) = 0D0
	afLocalForces(28) = 0D0
	afLocalForces(29) = 0D0
	afLocalForces(30) = 0D0
	afLocalForces(31) = 0D0
	afLocalForces(32) = 0D0
	afLocalForces(33) = 0D0
	afLocalForces(34) = 0D0
	afLocalForces(35) = 0D0
	afLocalForces(36) = 0D0
	afLocalForces(37) = 0D0
	afLocalForces(38) = 0D0
	afLocalForces(39) = 0D0
	afLocalForces(40) = 0D0
	afLocalForces(41) = 0D0
	afLocalForces(42) = 0D0
	afLocalForces(43) = 0D0
	afLocalForces(44) = 0D0
	afLocalForces(45) = 0D0
	afLocalForces(46) = 0D0
	afLocalForces(47) = 0D0
	afLocalForces(48) = 0D0
	afLocalForces(49) = 0D0
	afLocalForces(50) = 0D0
	afLocalForces(51) = 0D0
	afLocalForces(52) = 0D0
	afLocalForces(53) = 0D0
	afLocalForces(54) = 0D0
	afLocalForces(55) = 0D0
	afLocalForces(56) = 0D0
	afLocalForces(57) = 0D0
	afLocalForces(58) = 0D0
	afLocalForces(59) = 0D0
	afLocalForces(60) = 0D0
	If (bIsLumped) Then
		I = CalcH20MLumped(iInt, fDensity, afWeights, afM)
	Else
		I = CalcH20MConsistent(iInt, fDensity, afS, afWeights, afM)
	End If

! Multiply SKYLINE Element M with vector pfAcc and store to fLocalForces
	iDataPos = 1
	Do I = 1, 60
		iAccR = mod(I - 1, 3) + 1
		afLocalForces(I) = afLocalForces(I) + afM(iDataPos) * afAcc(iAccR)
		iDataPos = iDataPos + 1
		Do J = I + 1, 60
			afLocalForces(I) = afLocalForces(I) + afM(iDataPos) * afAcc(mod(J - 1, 3) + 1)
			afLocalForces(J) = afLocalForces(J) + afM(iDataPos) * afAcc(iAccR)
			iDataPos = iDataPos + 1
		End Do
	End Do
End Subroutine

Subroutine CalcH20ForcesFromAcc(piEqLM, pfElementM, pfAcc, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementM(1), pfAcc(3)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iRow, iAccR, iTemp, iDataPos
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

! Multiply SKYLINE Element M with vector pfAcc and store to fLocalForces
	iDataPos = 1
	Do I = 1, 60
		iAccR = mod(I - 1, 3) + 1
		fLocalForces(I) = fLocalForces(I) + pfElementM(iDataPos) * pfAcc(iAccR)
		iDataPos = iDataPos + 1
		Do J = I + 1, 60
			fLocalForces(I) = fLocalForces(I) + pfElementM(iDataPos) * pfAcc(mod(J - 1, 3) + 1)
			fLocalForces(J) = fLocalForces(J) + pfElementM(iDataPos) * pfAcc(iAccR)
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

