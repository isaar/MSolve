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
Subroutine CalcH8NablaShape(fXi, fEta, fZeta, faDS)
	Implicit None

	Real, Intent(IN) :: fXi, fEta, fZeta
	Real, Intent(OUT) :: faDS(24)
	Real, Parameter :: fSq125 = 0.35355339059327376220042218105242
	Real fXiP, fEtaP, fZetaP, fXiM, fEtaM, fZetaM

! Main routine
	fXiP = (1D0 + fXi) * fSq125
	fEtaP = (1D0 + fEta) * fSq125
	fZetaP = (1D0 + fZeta) * fSq125
	fXiM = (1D0 - fXi) * fSq125
	fEtaM = (1D0 - fEta) * fSq125
	fZetaM = (1D0 - fZeta) * fSq125

! Calculation of natural coordinate derivatives of the shape functions
! Corresponding to XI
	faDS(1) = - fEtaM * fZetaM
	faDS(2) = - faDS(1)
	faDS(3) = fEtaP * fZetaM
	faDS(4) = - faDS(3)
	faDS(5) = - fEtaM * fZetaP
	faDS(6) = - faDS(5)
	faDS(7) = fEtaP * fZetaP
	faDS(8) = - faDS(7)
! Corresponding to ETA
    faDS(9) = - fXiM * fZetaM
    faDS(10) = - fXiP * fZetaM
    faDS(11) = - faDS(10)
    faDS(12) = - faDS(9)
	faDS(13) = - fXiM * fZetaP
	faDS(14) = - fXiP * fZetaP
	faDS(15) = - faDS(14)
	faDS(16) = - faDS(13)
! Corresponding to ZETA
	faDS(17) = - fXiM * fEtaM
	faDS(18) = - fXiP * fEtaM
	faDS(19) = - fXiP * fEtaP
	faDS(20) = - fXiM * fEtaP
	faDS(21) = - faDS(17)
	faDS(22) = - faDS(18)
	faDS(23) = - faDS(19)
	faDS(24) = - faDS(20)
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
Subroutine CalcH8Shape(fXi, fEta, fZeta, faS)
	Implicit None

	Real, Intent(IN) :: fXi, fEta, fZeta
	Real, Intent(OUT) :: faS(8)
	Real, Parameter :: fSqC125 = 0.5D0
	Real fXiP, fEtaP, fZetaP, fXiM, fEtaM, fZetaM

! Main routine
	fXiP = (1.0 + fXi) * fSqC125
	fEtaP = (1.0 + fEta) * fSqC125
	fZetaP = (1.0 + fZeta) * fSqC125
	fXiM = (1.0 - fXi) * fSqC125
	fEtaM = (1.0 - fEta) * fSqC125
	fZetaM = (1.0 - fZeta) * fSqC125

	faS(1) = fXiM * fEtaM * fZetaM
	faS(2) = fXiP * fEtaM * fZetaM
	faS(3) = fXiP * fEtaP * fZetaM
	faS(4) = fXiM * fEtaP * fZetaM
	faS(5) = fXiM * fEtaM * fZetaP
	faS(6) = fXiP * fEtaM * fZetaP
	faS(7) = fXiP * fEtaP * fZetaP
	faS(8) = fXiM * fEtaP * fZetaP
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
Integer Function CalcH8JDetJ(faXYZ, faDS, faJ, fDetJ)
	Implicit None
	Include 'hCommon.h'

	Real, Intent(IN) :: faXYZ(3, 8), faDS(24)
	Real, Intent(OUT) :: faJ(3, 3), fDetJ
	Real fDet1, fDet2, fDet3

! Main routine
	CalcH8JDetJ = iErrNegDet
! Compute the Jacobian at (XI, ETA, ZETA)
	faJ(1, 1) = faDS(1) * faXYZ(1, 1) + &
		faDS(2) * faXYZ(1, 2) + faDS(3) * faXYZ(1, 3) + &
		faDS(4) * faXYZ(1, 4) + faDS(5) * faXYZ(1, 5) + &
		faDS(6) * faXYZ(1, 6) + faDS(7) * faXYZ(1, 7) + &
		faDS(8) * faXYZ(1, 8)
	faJ(1, 2) = faDS(1) * faXYZ(2, 1) + &
		faDS(2) * faXYZ(2, 2) + faDS(3) * faXYZ(2, 3) + &
		faDS(4) * faXYZ(2, 4) + faDS(5) * faXYZ(2, 5) + &
		faDS(6) * faXYZ(2, 6) + faDS(7) * faXYZ(2, 7) + &
		faDS(8) * faXYZ(2, 8)
	faJ(1, 3) = faDS(1) * faXYZ(3, 1) + &
		faDS(2) * faXYZ(3, 2) + faDS(3) * faXYZ(3, 3) + &
		faDS(4) * faXYZ(3, 4) + faDS(5) * faXYZ(3, 5) + &
		faDS(6) * faXYZ(3, 6) + faDS(7) * faXYZ(3, 7) + &
		faDS(8) * faXYZ(3, 8)
	faJ(2, 1) = faDS(9) * faXYZ(1, 1) + &
		faDS(10) * faXYZ(1, 2) + faDS(11) * faXYZ(1, 3) + &
		faDS(12) * faXYZ(1, 4) + faDS(13) * faXYZ(1, 5) + &
		faDS(14) * faXYZ(1, 6) + faDS(15) * faXYZ(1, 7) + &
		faDS(16) * faXYZ(1, 8)
	faJ(2, 2) = faDS(9) * faXYZ(2, 1) + &
		faDS(10) * faXYZ(2, 2) + faDS(11) * faXYZ(2, 3) + &
		faDS(12) * faXYZ(2, 4) + faDS(13) * faXYZ(2, 5) + &
		faDS(14) * faXYZ(2, 6) + faDS(15) * faXYZ(2, 7) + &
		faDS(16) * faXYZ(2, 8)
	faJ(2, 3) = faDS(9) * faXYZ(3, 1) + &
		faDS(10) * faXYZ(3, 2) + faDS(11) * faXYZ(3, 3) + &
		faDS(12) * faXYZ(3, 4) + faDS(13) * faXYZ(3, 5) + &
		faDS(14) * faXYZ(3, 6) + faDS(15) * faXYZ(3, 7) + &
		faDS(16) * faXYZ(3, 8)
	faJ(3, 1) = faDS(17) * faXYZ(1, 1) + &
		faDS(18) * faXYZ(1, 2) + faDS(19) * faXYZ(1, 3) + &
		faDS(20) * faXYZ(1, 4) + faDS(21) * faXYZ(1, 5) + &
		faDS(22) * faXYZ(1, 6) + faDS(23) * faXYZ(1, 7) + &
		faDS(24) * faXYZ(1, 8)
	faJ(3, 2) = faDS(17) * faXYZ(2, 1) + &
		faDS(18) * faXYZ(2, 2) + faDS(19) * faXYZ(2, 3) + &
		faDS(20) * faXYZ(2, 4) + faDS(21) * faXYZ(2, 5) + &
		faDS(22) * faXYZ(2, 6) + faDS(23) * faXYZ(2, 7) + &
		faDS(24) * faXYZ(2, 8)
	faJ(3, 3) = faDS(17) * faXYZ(3, 1) + &
		faDS(18) * faXYZ(3, 2) + faDS(19) * faXYZ(3, 3) + &
		faDS(20) * faXYZ(3, 4) + faDS(21) * faXYZ(3, 5) + &
		faDS(22) * faXYZ(3, 6) + faDS(23) * faXYZ(3, 7) + &
		faDS(24) * faXYZ(3, 8)

	fDet1 = faJ(1, 1) * (faJ(2, 2) * faJ(3, 3) - faJ(3, 2) * faJ(2, 3))
	fDet2 = -faJ(1, 2) * (faJ(2, 1) * faJ(3, 3) - faJ(3, 1) * faJ(2, 3))
	fDet3 = faJ(1, 3) * (faJ(2, 1) * faJ(3, 2) - faJ(3, 1) * faJ(2, 2))
	fDetJ = fDet1 + fDet2 + fDet3
! Error if determinant is negative (return)
	If (fDetJ.LT.0.00000001) then
		Write(*, *) "Negative determinant (H8)"
		Return
	EndIf
	CalcH8JDetJ = iErrNoError
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
Subroutine CalcH8B(fDetJ, faJ, faDS, faB)
	Implicit None

	Real, Intent(IN) :: fDetJ, faJ(3, 3), faDS(24)
	Real, Intent(OUT) :: faB(6, 24)
	Real faJInv(3, 3), fDetInv
!
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
	faB(1, 1) = faJInv(1, 1) * faDS(1) + faJInv(1, 2) * faDS(9) + faJInv(1, 3) * faDS(17)
	faB(1, 2) = 0D0
	faB(1, 3) = 0D0
	faB(1, 4) = faJInv(1, 1) * faDS(2) + faJInv(1, 2) * faDS(10) + faJInv(1, 3) * faDS(18)
	faB(1, 5) = 0D0
	faB(1, 6) = 0D0
	faB(1, 7) = faJInv(1, 1) * faDS(3) + faJInv(1, 2) * faDS(11) + faJInv(1, 3) * faDS(19)
	faB(1, 8) = 0D0
	faB(1, 9) = 0D0
	faB(1, 10) = faJInv(1, 1) * faDS(4) + faJInv(1, 2) * faDS(12) + faJInv(1, 3) * faDS(20)
	faB(1, 11) = 0D0
	faB(1, 12) = 0D0
	faB(1, 13) = faJInv(1, 1) * faDS(5) + faJInv(1, 2) * faDS(13) + faJInv(1, 3) * faDS(21)
	faB(1, 14) = 0D0
	faB(1, 15) = 0D0
	faB(1, 16) = faJInv(1, 1) * faDS(6) + faJInv(1, 2) * faDS(14) + faJInv(1, 3) * faDS(22)
	faB(1, 17) = 0D0
	faB(1, 18) = 0D0
	faB(1, 19) = faJInv(1, 1) * faDS(7) + faJInv(1, 2) * faDS(15) + faJInv(1, 3) * faDS(23)
	faB(1, 20) = 0D0
	faB(1, 21) = 0D0
	faB(1, 22) = faJInv(1, 1) * faDS(8) + faJInv(1, 2) * faDS(16) + faJInv(1, 3) * faDS(24)
	faB(1, 23) = 0D0
	faB(1, 24) = 0D0
	faB(2, 1) = 0D0
	faB(2, 2) = faJInv(2, 1) * faDS(1) + faJInv(2, 2) * faDS(9) + faJInv(2, 3) * faDS(17)
	faB(2, 3) = 0D0
	faB(2, 4) = 0D0
	faB(2, 5) = faJInv(2, 1) * faDS(2) + faJInv(2, 2) * faDS(10) + faJInv(2, 3) * faDS(18)
	faB(2, 6) = 0D0
	faB(2, 7) = 0D0
	faB(2, 8) = faJInv(2, 1) * faDS(3) + faJInv(2, 2) * faDS(11) + faJInv(2, 3) * faDS(19)
	faB(2, 9) = 0D0
	faB(2, 10) = 0D0
	faB(2, 11) = faJInv(2, 1) * faDS(4) + faJInv(2, 2) * faDS(12) + faJInv(2, 3) * faDS(20)
	faB(2, 12) = 0D0
	faB(2, 13) = 0D0
	faB(2, 14) = faJInv(2, 1) * faDS(5) + faJInv(2, 2) * faDS(13) + faJInv(2, 3) * faDS(21)
	faB(2, 15) = 0D0
	faB(2, 16) = 0D0
	faB(2, 17) = faJInv(2, 1) * faDS(6) + faJInv(2, 2) * faDS(14) + faJInv(2, 3) * faDS(22)
	faB(2, 18) = 0D0
	faB(2, 19) = 0D0
	faB(2, 20) = faJInv(2, 1) * faDS(7) + faJInv(2, 2) * faDS(15) + faJInv(2, 3) * faDS(23)
	faB(2, 21) = 0D0
	faB(2, 22) = 0D0
	faB(2, 23) = faJInv(2, 1) * faDS(8) + faJInv(2, 2) * faDS(16) + faJInv(2, 3) * faDS(24)
	faB(2, 24) = 0D0
	faB(3, 1) = 0D0
	faB(3, 2) = 0D0
	faB(3, 3) = faJInv(3, 1) * faDS(1) + faJInv(3, 2) * faDS(9) + faJInv(3, 3) * faDS(17)
	faB(3, 4) = 0D0
	faB(3, 5) = 0D0
	faB(3, 6) = faJInv(3, 1) * faDS(2) + faJInv(3, 2) * faDS(10) + faJInv(3, 3) * faDS(18)
	faB(3, 7) = 0D0
	faB(3, 8) = 0D0
	faB(3, 9) = faJInv(3, 1) * faDS(3) + faJInv(3, 2) * faDS(11) + faJInv(3, 3) * faDS(19)
	faB(3, 10) = 0D0
	faB(3, 11) = 0D0
	faB(3, 12) = faJInv(3, 1) * faDS(4) + faJInv(3, 2) * faDS(12) + faJInv(3, 3) * faDS(20)
	faB(3, 13) = 0D0
	faB(3, 14) = 0D0
	faB(3, 15) = faJInv(3, 1) * faDS(5) + faJInv(3, 2) * faDS(13) + faJInv(3, 3) * faDS(21)
	faB(3, 16) = 0D0
	faB(3, 17) = 0D0
	faB(3, 18) = faJInv(3, 1) * faDS(6) + faJInv(3, 2) * faDS(14) + faJInv(3, 3) * faDS(22)
	faB(3, 19) = 0D0
	faB(3, 20) = 0D0
	faB(3, 21) = faJInv(3, 1) * faDS(7) + faJInv(3, 2) * faDS(15) + faJInv(3, 3) * faDS(23)
	faB(3, 22) = 0D0
	faB(3, 23) = 0D0
	faB(3, 24) = faJInv(3, 1) * faDS(8) + faJInv(3, 2) * faDS(16) + faJInv(3, 3) * faDS(24)
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
Subroutine CalcH8GaussMatrices(iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8GaussMatrices
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faXYZ(3, 8)
	Real, Intent(OUT) :: faDS(24, iInt ** 3), faS(8, iInt ** 3), &
		faB(6, 24, iInt ** 3), faDetJ(iInt ** 3), faJ(3, 3, iInt ** 3), &
		faWeight(iInt ** 3)
	Real fXi, fEta, fZeta
	Integer M, iX, iY, iZ, iErr, CalcH8JDetJ
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
				Call CalcH8NablaShape(fXi, fEta, fZeta, faDS(1, M))
				Call CalcH8Shape(fXi, fEta, fZeta, faS(1, M))
				iErr = CalcH8JDetJ(faXYZ, faDS(1, M), faJ(1, 1, M), faDetJ(M))
				Call CalcH8B(faDetJ(M), faJ(1, 1, M), faDS(1, M), faB(1, 1, M))
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
Subroutine CalcH8K(iInt, faE, faB, faWeight, faK)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8K
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faWeight(iInt ** 3), faB(6, 24, iInt ** 3), faE(6, 6, iInt ** 3)
	Real, Intent(OUT) :: faK(300)
	Real faEB(6)
	Real fStiff, fE1, fE2, fE3, fE4
	Integer I, J, K, M, KS, iX, iY, iZ

! Main routine 
! Calculate stiffness matrix
	M = 0
	Do I = 1, 300
		faK(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
! Add contribution to element stiffness
				KS = 0
				Do I = 1, 24
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

					Do J = I, 24
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
Integer Function CalcH8MConsistent(iInt, fDensity, faS, faWeight, faM)
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fDensity, faWeight(iInt ** 3), faS(8, iInt ** 3)
	Real, Intent(OUT) :: faM(300)
	Integer I, J, K, M, KS, iX, iY, iZ, iWidth, iS2, iS1
	Real fWeight

! Main routine
	M = 0
	Do I = 1, 300
		faM(I) = 0D0
	EndDo
	Do iX = 1, iInt
		Do iY = 1, iInt
			Do iZ = 1, iInt
				M = M + 1
				KS = -2
				fWeight = faWeight(M) * fDensity
				Do I = 1, 8
					Do J = I, 8
						KS = KS + 3
						faM(KS) = faM(KS) + faS(I, M) * faS(J, M) * fWeight
					EndDo

					iWidth = 27 - 3 * I
					Do J = I, 8
						KS = KS + 3
						faM(KS) = faM(KS - iWidth)
					EndDo
					KS = KS - 1
					iWidth = 26 - 3 * I
					Do J = I, 8
						KS = KS + 3
						faM(KS) = faM(KS - iWidth)
					EndDo
					KS = KS - 2
				EndDo
			EndDo
		EndDo
	EndDo
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
Subroutine CalcH8Strains(iInt, faB, fau, faStrains)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8Strains
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faB(6, 24, iInt ** 3), fau(24)
	Real, Intent(OUT) :: faStrains(6, iInt ** 3)
	Integer :: M

! Main routine 
	Do M = 1, iInt ** 3
		faStrains(:, M) = MatMul(faB(:, :, M), fau)
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
Integer Function CalcH8MLumped(iInt, fDensity, faWeight, faM)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8MLumped
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: fDensity, faWeight(iInt ** 3)
	Real, Intent(OUT) :: faM(300)
	Integer M, iS2, iS1
	Real fMassTerm

! Main routine
	Do M = 1, 300
		faM(M) = 0D0
	EndDo

	fMassTerm = 0D0
	Do M = 1, iInt ** 3
		fMassTerm = fMassTerm + faWeight(M)
	EndDo	
	fMassTerm = fMassTerm * fDensity * 0.125D0

	iS2 = 1
	Do iS1 = 24, 1, -1
		faM(iS2) = fMassTerm
		iS2 = iS2 + iS1
	EndDo
End Function

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
Subroutine CalcH8Forces(iInt, faB, faWeight, faStresses, faForces)
  !DEC$ ATTRIBUTES DLLEXPORT::CalcH8Forces
	Implicit None

	Integer, Intent(IN) :: iInt
	Real, Intent(IN) :: faWeight(iInt ** 3), faB(6, 24, iInt ** 3), faStresses(6, iInt ** 3)
	Real, Intent(OUT) :: faForces(24)
	Integer :: M, iX, iY, iZ

! Main routine 
! Calculate stiffness matrix
	M = 0
	Do iX = 1, 24
		faForces(iX) = 0D0
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

Subroutine CalcH8ForcesFromDisp(piEqLM, pfElementK, pfU, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementK(1), pfU(24)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iRow, iAccR, iTemp, iDataPos
	Real :: fLocalForces(24)

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

! Multiply SKYLINE Element K with vector pfU and store to fLocalForces
	iDataPos = 1
	Do I = 1, 24
		fLocalForces(I) = fLocalForces(I) + pfElementK(iDataPos) * pfU(I)
		iDataPos = iDataPos + 1
		Do J = I + 1, 24
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
End Subroutine

Subroutine CalcH8ForcesAcc(iInt, bIsLumped, fDensity, afAcc, afS, afWeights, afLocalForces)
	Implicit None
	Integer, Intent(IN) :: iInt
	Logical, Intent(IN) :: bIsLumped
	Real, Intent(IN) :: fDensity, afS(8, iInt ** 3), afAcc(3), afWeights(*)
	Real, Intent(OUT) :: afLocalForces(24)
	Integer :: I, J, iAccR, iDataPos, CalcH8MLumped, CalcH8MConsistent
	Real :: afM(300)

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
	If (bIsLumped) Then
		I = CalcH8MLumped(iInt, fDensity, afWeights, afM)
	Else
		I = CalcH8MConsistent(iInt, fDensity, afS, afWeights, afM)
	End If

! Multiply SKYLINE Element M with vector pfAcc and store to fLocalForces
	iDataPos = 1
	Do I = 1, 24
		iAccR = mod(I - 1, 3) + 1
		afLocalForces(I) = afLocalForces(I) + afM(iDataPos) * afAcc(iAccR)
		iDataPos = iDataPos + 1
		Do J = I + 1, 24
			afLocalForces(I) = afLocalForces(I) + afM(iDataPos) * afAcc(mod(J - 1, 3) + 1)
			afLocalForces(J) = afLocalForces(J) + afM(iDataPos) * afAcc(iAccR)
			iDataPos = iDataPos + 1
		End Do
	End Do
End Subroutine

Subroutine CalcH8ForcesFromAcc(piEqLM, pfElementM, pfAcc, pfGlobalForces)
	Implicit None

	Integer, Intent(IN) :: piEqLM(1)
	Real, Intent(IN) :: pfElementM(1), pfAcc(3)
	Real, Intent(OUT) :: pfGlobalForces(1)

	Integer :: I, J, iRow, iAccR, iTemp, iDataPos
	Real :: fLocalForces(24)

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

! Multiply SKYLINE Element M with vector pfAcc and store to fLocalForces
	iDataPos = 1
	Do I = 1, 24
		iAccR = mod(I - 1, 3) + 1
		fLocalForces(I) = fLocalForces(I) + pfElementM(iDataPos) * pfAcc(iAccR)
		iDataPos = iDataPos + 1
		Do J = I + 1, 24
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
End Subroutine

