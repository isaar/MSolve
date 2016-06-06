! ******************************************
! File: hCommon.for
!
! Purpose: Contains common parameters
! ******************************************

! Sides
	Integer, Parameter :: iXZDown = 1
	Integer, Parameter :: iXZUp = 2
	Integer, Parameter :: iYZLeft = 3
	Integer, Parameter :: iYZRight = 4
	Integer, Parameter :: iXYFront = 5
	Integer, Parameter :: iXYBack = 6

! Error log file ID
	Parameter iOutError = 20

! Error codes
	Parameter iErrNoError = 0
	Parameter iErrUnknown = -1
	Parameter iErrFileNotFound = -10
	Parameter iErrReadFormat = -11
	Parameter iErrNameExists = -12
	Parameter iErrNodeIxNoMatch = -50
	Parameter iErrNameNoMatch = -51
	Parameter iErrNotEnoughMemory = -100
	Parameter iErrZeroNodes = -200
	Parameter iErrUnknownElementType = -210
	Parameter iErrNegDet = -1001

! Group codes
	Parameter iGrpFirst = 1
	Parameter iGrpModelData = 1
	Parameter iGrpElementData = 2
	Parameter iGrpLoadData = 3
	Parameter iGrpSolverData = 4
	Parameter iGrpLast = 4

! Parser codes
	Parameter iComRemarks = 99				!** Any kind of remarks
	Parameter iComAnalysisStatic = 1		!*Analysis, Type=Static
	Parameter iComAnalysisDynNewmark = 2	!*Analysis, Type=Dynamic NEWMARK, Mass=xxx
	Parameter iComAnalysisDynExplicit = 3	!*Analysis, Type=Dynamic Explicit, Mass=xxx
	Parameter iComNodes = 10				!*Nodes
	Parameter iComConstraints = 20			!*Constraints
	Parameter iComProperty = 30				!*Property, Type=1, Name=xxx
	Parameter iComElements = 40				!*Elements, Type=x, Props=xxx
	Parameter iComCurve = 50				!*Curve, name=xxx
	Parameter iComLoadStaticNodal = 61		!*Loads, Type=Static Nodal, Name=xxx
	Parameter iComLoadDynamicNodal = 62		!*Loads, Type=Dynamic Nodal curve, Name=xxx
	Parameter iComLoadDynamicNodalFile = 63	!*Loads, Type=Dynamic file, Name=xxx
	Parameter iComLoadDynamicNodalAcc = 64	!*Loads, Type=Dynamic file accelogram, Name=xxx
	Parameter iComLoadStaticGlobalAcc = 65  !*Loads, Type=Global acceleration vector, Name=xxx
	Parameter iComLoadStaticNodalDisp = 66  !*Loads, Type=Static Nodal DOF Value, Name=xxx
	Parameter iComLoadCase = 70				!*Load case, Name=xxx
	Parameter iComInitialState = 71			!*Initial state
	Parameter iComOutput = 80				!*Output
	Parameter iComRayleigh = 90				!*Rayleigh
	Parameter iComSubdomain = 100			!*Subdomain
	Parameter iComProcessor = 101			!*Processor

! Property types
	Parameter iPropPlate1 = 16
	Parameter iPropPlate2 = 17
	Parameter iPropPlaneStrain1 = 18
	Parameter iPropPlaneStrain2 = 19

! Element codes
	Parameter iElmHexa8 = 2
	Parameter iElmHexa20 = 4
	Parameter iElmCST = 11
	Parameter iElmWedge15 = 51
	Parameter iElmHexa8u8p = 101
	Parameter iElmHexa20u8p = 102

! DOF codes
	Parameter iDOFuX = 1
	Parameter iDOFuY = 2
	Parameter iDOFuZ = 3
	Parameter iDOFp = 11
!	Parameter iDOFpX = 11
!	Parameter iDOFpY = 12
!	Parameter iDOFpZ = 13
	Parameter iSpringX = 21
	Parameter iSpringY = 22
	Parameter iSpringZ = 23
	Parameter iDOFFluid = 30

! Solver codes
	Parameter iSolverSkyline = 1
	Parameter iSolverSkylineBB = 2
	Parameter iSolverFETI_T = 10
	Parameter iSolverPCG = 20
	Parameter iSolverPCG_Matlab = 101
	Parameter iSolverPSM_Matlab = 102
	Parameter iSolverFETI_Matlab = 110
	Parameter iSolverPCG_Porous_Matlab = 121
	Parameter iSolverFETI_NoCoarse = 200

! FETI solver codes
	Parameter iFSFETI = 1
	Parameter iFSFETI_T = 10

! Preconditioner codes
	Parameter iPrecNone = 0
	Parameter iPrecLumped = 1
	Parameter iPrecDirichlet = 2
	Parameter iPrecDiagonal = 11

! Property codes
	Parameter iPrpPorous = 101

! Material codes
	Parameter iMatElastic = 1
	Parameter iMatVonMises = 2
	Parameter iMatVonMisesConsistent = 3
	Parameter iMatIsotropicDamageLinear = 10
	Parameter iMatIsotropicDamageExp = 11
	Parameter iMatTenProperties = 9

! Non-linear solvers
	Parameter iNLSNewton = 1
	Parameter iNLSNewtonModified = 2

! Non-linear output strategy
	Parameter iNLODefault = 0
	Parameter iNLOIntermediate = 1
	Parameter iNLOConverged = 2
	Parameter iNLOAll = 3

! Matrix codes
	Parameter iMatrixK = 1
	Parameter iMatrixM = 2
	Parameter iMatrixC = 3
	Parameter iMatrixH = 4
	Parameter iMatrixS = 5
	Parameter iMatrixQ = 6
	Parameter iMatrixQTilde = 7
	Parameter iMatrixShape = 11
	Parameter iMatrixB = 12
