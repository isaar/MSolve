using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Primary multiscale analysis class that connects all nesessary structures for a FE2 simulation for 3D continuum structures
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class MicrostructureDefGrad3D : StructuralProblemsMicrostructureBase_v2, IContinuumMaterial3DDefGrad_v2
    {
        private Model_v2 model { get; set; }
        //private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private Dictionary<int, Node_v2> boundaryNodes { get; set; }
        Dictionary<int, Dictionary<int, Element_v2>> boundaryElements;
        private IRVEbuilder_v2 rveBuilder;
        private bool EstimateOnlyLinearResponse;
        //private NewtonRaphsonNonLinearAnalyzer microAnalyzer;
        private double volume;
        public Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain { get; private set; }
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        private IScaleTransitions_v2 scaleTransitions = new DefGradVec3DScaleTransition_v2();
        Random rnd1 = new Random();
        private readonly Func<Model_v2, ISolver_v2> createSolver;

        // aparaithta gia to implementation tou IFiniteElementMaterial3D
        Matrix constitutiveMatrix;
        private double[] SPK_vec=new double[6];
        private bool modified; // opws sto MohrCoulomb gia to modified

        private double[,] Cijrs_prev;
        private bool matrices_not_initialized = true;
        private double tol;
        public void InitializeMatrices()
        {
            Cijrs_prev = new double[6, 6];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
            constitutiveMatrix = Matrix.CreateZero(6, 6);
        }


        //double[] Stresses { get; }
        //IMatrix2D ConstitutiveMatrix { get; } TODOGerasimos

        //Random properties 
        private int database_size;

        public MicrostructureDefGrad3D(IRVEbuilder_v2 rveBuilder, Func<Model_v2, ISolver_v2> createSolver, 
            bool EstimateOnlyLinearResponse, int database_size)
        {
            this.rveBuilder = rveBuilder;
            this.createSolver = createSolver;
            this.EstimateOnlyLinearResponse = EstimateOnlyLinearResponse;
            this.database_size = database_size;            
        }

        private void InitializeData()
        {
            Tuple<Model_v2, Dictionary<int, Node_v2>, double> modelAndBoundaryNodes = this.rveBuilder.GetModelAndBoundaryNodes();
            this.model = modelAndBoundaryNodes.Item1;
            this.boundaryNodes = modelAndBoundaryNodes.Item2;
            this.boundaryElements = GetSubdomainsBoundaryFiniteElementsDictionaries_v2(model, boundaryNodes);
            this.volume = modelAndBoundaryNodes.Item3;
            DefineAppropriateConstraintsForBoundaryNodes();
            this.model.ConnectDataStructures();            
        }


        private void DefineAppropriateConstraintsForBoundaryNodes()
        {
            foreach(Node_v2 boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        private void InitializeFreeAndPrescribedDofsInitialDisplacementVectors()
        {
            uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, IVector>();
            foreach(Subdomain_v2 subdomain in model.SubdomainsDictionary.Values)
            {
                uInitialFreeDOFDisplacementsPerSubdomain.Add(subdomain.ID, Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs));// prosoxh sto Id twn subdomain
            }
            double[,] DGtr = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] DefGradVec = new double[9] { DGtr[0, 0], DGtr[1, 1], DGtr[2, 2], DGtr[1, 0], DGtr[2, 1], DGtr[0, 2], DGtr[2, 0], DGtr[0, 1], DGtr[1, 2], };
            initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                DefGradVec, initialConvergedBoundaryDisplacements);
            }            
        }

        public object Clone()
        {
            int new_rve_id = rnd1.Next(1, database_size + 1);
            return new MicrostructureDefGrad3D(rveBuilder.Clone(new_rve_id), createSolver, EstimateOnlyLinearResponse, database_size);
        }

        public Dictionary<int, Node_v2> BoundaryNodesDictionary
        {
            get { return boundaryNodes; }
        }
        public IList<Node_v2> BoundaryNodes
        {
            get { return boundaryNodes.Values.ToList<Node_v2>(); }
        }

        public void UpdateMaterial(double[] DefGradVec)
        {
            ISolver_v2 solver;
            if (matrices_not_initialized)
            {
                this.InitializeMatrices();
                this.InitializeData();
                solver = createSolver(model);
                solver.OrderDofs(false);
                foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
                {
                    linearSystem.Reset(); //TODO find out if new structures cause any problems
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
                this.InitializeFreeAndPrescribedDofsInitialDisplacementVectors();
            }
            else
            {
                solver = createSolver(model);
                solver.OrderDofs(false); //v2.1. TODO: Is this needed in this case?
                foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
                {
                    linearSystem.Reset();
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces; //TODO MS 
                }
            }

            for (int i1 = 0; i1 < 6; i1++)
            {
                for (int j1 = 0; j1 < 6; j1++)
                {Cijrs_prev[i1, j1] = constitutiveMatrix[i1, j1];}
            }

            #region Rve prescribed Dofs total DIsplacement Dictionary Creation (nessesary for NRNLAnalyzer)
            Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                DefGradVec, totalPrescribedBoundaryDisplacements);
            }
            #endregion
                     
            //var linearSystems = CreateNecessaryLinearSystems(model);    // OPOU pairnei rhs apo subdomainForces       
            //var solver = GetAppropriateSolver(linearSystems);

            
            #region Creation of nessesary analyzers for NRNLAnalyzer and Creation of Microstructure analyzer (NRNLdevelop temporarilly) and solution ;
            int increments = 1; int MaxIterations = 100; int IterationsForMatrixRebuild = 1;
            (MicrostructureBvpNRNLAnalyzer microAnalyzer, ProblemStructural_v2 provider, ElementStructuralStiffnessProvider_v2 elementProvider) = 
                AnalyzeMicrostructure_v2(model, solver, increments, MaxIterations, IterationsForMatrixRebuild,
                totalPrescribedBoundaryDisplacements, initialConvergedBoundaryDisplacements, boundaryNodes, uInitialFreeDOFDisplacementsPerSubdomain);
            #endregion

            #region update of free converged displacements vectors
            uInitialFreeDOFDisplacementsPerSubdomain = microAnalyzer.GetConvergedSolutionVectorsOfFreeDofs();// ousiastika to u pou twra taftizetai me to uPlusuu
            #endregion


            #region INTEGRATION stresses 
            Dictionary<int, IVector> du = microAnalyzer.GetConvergedIncrementalSolutionVectorsOfFreeDofs();
            Dictionary<int, double[]> FppReactionVectorSubdomains = SubdomainCalculationsMultiple_v2.CalculateFppReactionsVectorSubdomains_v2(model, elementProvider, scaleTransitions, boundaryNodes,
                uInitialFreeDOFDisplacementsPerSubdomain, du, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, increments, increments);
            double[] FppReactionVector= SubdomainCalculationsMultiple_v2.CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal_v2(FppReactionVectorSubdomains);



            double[] DqFpp = SubdomainCalculations_v2.CalculateDqFpp_v2(FppReactionVector, scaleTransitions, boundaryNodes);

            double[] FPK_vec = new double [DqFpp.Length];
            for (int i1 = 0; i1 < DqFpp.Length; i1++)
            { FPK_vec[i1]=(1 / volume) * DqFpp[i1]; }

            double[,] DefGradMat = new double[3, 3] { { DefGradVec[0], DefGradVec[3], DefGradVec[6] }, { DefGradVec[7], DefGradVec[1], DefGradVec[4] }, { DefGradVec[5], DefGradVec[8], DefGradVec[2] } };
            double[,] FPK_mat = new double[3, 3] { { FPK_vec[0], FPK_vec[3], FPK_vec[6] }, { FPK_vec[7], FPK_vec[1], FPK_vec[4] }, { FPK_vec[5], FPK_vec[8], FPK_vec[2] } };
            double[,] SPK_mat = transformFPKtoSPK(DefGradMat, FPK_mat);
            SPK_vec = new double[6] { SPK_mat[0,0], SPK_mat[1,1], SPK_mat[2,2], SPK_mat[0,1], SPK_mat[1,2], SPK_mat[0,2] };
            //TODOna elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            #endregion

            #region INTEGRATION constitutive Matrix
            var integrationSimultaneous = new SubdomainCalculationsAndAssembly();
            (Dictionary<int, double[][]> KfpDqSubdomains, Dictionary<int, double[][]> KppDqVectorsSubdomains) = 
                integrationSimultaneous.UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje_v2(model, elementProvider, scaleTransitions, boundaryNodes, boundaryElements, solver);


            Dictionary<int, double[][]> f2_vectorsSubdomains = SubdomainCalculationsMultiple_v2.CalculateKffinverseKfpDqSubdomains_v2(KfpDqSubdomains, model, elementProvider, scaleTransitions, boundaryNodes, solver);

            Dictionary<int, double[][]> f3_vectorsSubdomains = SubdomainCalculationsMultiple_v2.CalculateKpfKffinverseKfpDqSubdomains_v2(f2_vectorsSubdomains, model, elementProvider, scaleTransitions, boundaryNodes);

            double[][] f3_vectors = SubdomainCalculationsMultiple_v2.CombineMultipleSubdomainsIntegrationVectorsIntoTotal_v2(f3_vectorsSubdomains,scaleTransitions);
            double[][] KppDqVectors = SubdomainCalculationsMultiple_v2.CombineMultipleSubdomainsIntegrationVectorsIntoTotal_v2(KppDqVectorsSubdomains,scaleTransitions);

            double[][] f4_vectors = SubdomainCalculations_v2.SubtractConsecutiveVectors_v2(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations_v2.CalculateDqCondDq_v2(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] d2W_dfdf = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    d2W_dfdf[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            #endregion

            #region update of prescribed converged displacements vectors;
            initialConvergedBoundaryDisplacements = totalPrescribedBoundaryDisplacements;
            #endregion

            #region constitutive tensors transformation methods
            double[,] d2W_dFtrdFtr = Reorder_d2Wdfdf_to_d2W_dFtrdFtr(d2W_dfdf);

            double[,] Cinpk = Transform_d2WdFtrdFtr_to_Cijrs(d2W_dFtrdFtr, SPK_mat, DefGradMat); // to onomazoume Cinpk epeidh einai to 9x9 kai to diakrinoume etsi apo to Cijrs 6x6
            
            double[,] Cijrs = CombineCinpkTensorTermsIntoMatrix(Cinpk);
            
            #endregion

            constitutiveMatrix = Matrix.CreateFromArray(Cijrs);

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged(); 
        }        

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (Math.Abs(Cijrs_prev[i, j] - constitutiveMatrix[i, j]) > 1e-10)
                        return true;

            return false;
        }



        #region IFiniteElementMaterial3D methodoi mia mia 

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) CalculateOriginalConstitutiveMatrixWithoutNLAnalysis(); // TODOGerasimos arxiko constitutive mporei na upologizetai pio efkola
                return constitutiveMatrix;
            }
        }

        public double[] Stresses // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return SPK_vec; }
        }

        public void SaveState()
        {
            var subdomainUpdaters = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions_v2>(1); 
            foreach (Subdomain_v2 subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            {
                subdomainUpdaters.Add(subdomain.ID, new NonLinearSubdomainUpdaterWithInitialConditions_v2(subdomain)); //v2.3
                //subdomainUpdaters[counter] = new NonLinearSubdomainUpdaterWithInitialConditions(subdomain);
            }
            foreach (var subdomainUpdater in subdomainUpdaters.Values)
            {
                subdomainUpdater.UpdateState();
            }
            
        }

        public bool Modified
        {
            get { return modified; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public int ID
        {
            get { return 1000; }
        }


        #endregion
        // methodoi ews edw xrhsimopoiountai
        public void ClearState() 
        {
            // pithanws TODO 
        }
        public void ClearStresses()
        {
            // pithanws TODO 
        }
        public double[] Coordinates { get; set; }

        public double YoungModulus => throw new NotSupportedException();

        public double PoissonRatio => throw new NotSupportedException();


        #region transformation methods
        private double[,] transformFPKtoSPK(double[,] DefGradMat, double[,] FPK_mat)
        {
            double[,] SPK_Mat = new double[3, 3];
            //Vector solution = new Vector(new double[3]);

            for (int j1 = 0; j1 < 3; j1++)
            {
                double[] RHS = new double[3] { FPK_mat[0, j1], FPK_mat[1, j1], FPK_mat[2, j1] };
                double[] solution = commonCalculations.invert3by3(DefGradMat, RHS);                
                for (int i1 = 0; i1 < 3; i1++)
                {
                    SPK_Mat[i1, j1] = solution[i1];
                }
            }

            return SPK_Mat;
        }

        
        private double[,] Transform_d2WdFtrdFtr_to_Cijrs(double[,] Aijkl, double[,] SPK, double[,] F)
        {
            int[,] i_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };
            int[,] k_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };

            double[,] Cinpk = new double[9, 9];

            double[,] F__F__ = new double[9, 9];
            

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            F__F__[3 * i1 + k, 3 * j1 + l] = F[k, j1] * F[i1, l];
                        }
                    }

                }
            }

            Matrix F__F__Mat = Matrix.CreateFromArray(F__F__);            

            double[,] multipleRHSs=new double [9,9];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    double[] A_j_l = new double[9] { Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1]};

                    double[] sec_term = new double[9] { -SPK[i1 , k1 ], 0, 0, 0, -SPK[i1 , k1 ], 0, 0, 0, -SPK[i1, k1] };
                    
                    int RHScolumn = i1 * 3 + k1;
                    for (int a1 = 0; a1 < 9; a1++)
                    {
                        multipleRHSs[a1, RHScolumn] = A_j_l[a1] + sec_term[a1];
                    }
                    
                }
            }

            //TODO use solution multiple RHSs when pavailable :Matrix2D MultipleSolutions = F__F__Mat.SolveLU(new Matrix2D(multipleRHSs), true);
            Matrix inverse = F__F__Mat.Invert();
            Matrix MultipleSolutions = Matrix.CreateZero(9, 9);
            for (int i1=0; i1<9; i1++)
            {
                double[] RHS = new double[9];
                for (int i2 = 0; i2 < 9; i2++)
                {
                    RHS[i2] = multipleRHSs[i2, i1];
                }
                Vector solution = inverse * Vector.CreateFromArray(RHS);
                for (int i2 = 0; i2 < 9; i2++)
                {
                    MultipleSolutions[i2,i1] = solution[i2];
                }

            }


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {                    
                    int RHScolumn = i1 * 3 + k1;                    

                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[0,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[1,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[2,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[3,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[4,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[5,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[6,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[7,RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[8,RHScolumn];

                }
            }


            return Cinpk;
        }

        private double[,] CombineCinpkTensorTermsIntoMatrix(double[,] Cinpk)
        {
            // transformation se 6x6 se 2 vhmata

            double[,] Cijrs_columns = new double[9, 6];
            for (int i1 = 0; i1 < 9; i1++)
            {
                Cijrs_columns[i1, 0] = Cinpk[i1, 0];
                Cijrs_columns[i1, 1] = Cinpk[i1, 1];
                Cijrs_columns[i1, 2] = Cinpk[i1, 2];
                Cijrs_columns[i1, 3] = 0.5 * (Cinpk[i1, 3] + Cinpk[i1, 7]);
                Cijrs_columns[i1, 4] = 0.5 * (Cinpk[i1, 4] + Cinpk[i1, 8]);
                Cijrs_columns[i1, 5] = 0.5 * (Cinpk[i1, 5] + Cinpk[i1, 6]);
            }

            double[,] Cijrs = new double[6, 6];

            for (int j1 = 0; j1 < 6; j1++)
            {
                Cijrs[0, j1] = Cijrs_columns[0, j1];
                Cijrs[1, j1] = Cijrs_columns[1, j1];
                Cijrs[2, j1] = Cijrs_columns[2, j1];
                Cijrs[3, j1] = 0.5 * (Cijrs_columns[3, j1] + Cijrs_columns[7, j1]);
                Cijrs[4, j1] = 0.5 * (Cijrs_columns[4, j1] + Cijrs_columns[8, j1]);
                Cijrs[5, j1] = 0.5 * (Cijrs_columns[5, j1] + Cijrs_columns[6, j1]);
            }

            return Cijrs;
        }

        private double[,] Reorder_d2Wdfdf_to_d2W_dFtrdFtr(double[,] d2W_dfdf)
        {
            int[,] matLineData = new int[3, 3] { { 1, 4, 7 }, { 8, 2, 5 }, { 6, 9, 3 } };

            double[,] d2W_dFtrdFtr = new double[9, 9];

            for (int i1 = 1; i1 < 4; i1++)
            {
                for (int i2 = 1; i2 < 4; i2++)
                {
                    for (int i3 = 1; i3 < 4; i3++)
                    {
                        for (int i4 = 1; i4 < 4; i4++)
                        {

                            int d2 = i1; int d1 = i2; int d4 = i3; int d3 = i4;

                            int matLineA = matLineData[i1 - 1, i2 - 1]; //meion 1 logw zero based
                            int matRowA = matLineData[i3 - 1, i4 - 1];
                            int matLineW = matLineData[d1 - 1, d2 - 1];
                            int matRowW = matLineData[d3 - 1, d4 - 1];

                            d2W_dFtrdFtr[matLineA - 1, matRowA - 1] = d2W_dfdf[matLineW - 1, matRowW - 1];
                        }
                    }
                }
            }
            return d2W_dFtrdFtr;
        }
        #endregion

        public void CalculateOriginalConstitutiveMatrixWithoutNLAnalysis()
        {
            ISolver_v2 solver;
            if (matrices_not_initialized)
            {
                this.InitializeMatrices();
                this.InitializeData();
                solver = createSolver(model);
                solver.OrderDofs(false);
                foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
                {
                    linearSystem.Reset(); //TODO find out if new structures cause any problems
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
                this.InitializeFreeAndPrescribedDofsInitialDisplacementVectors();
            }
            else
            {
                solver = createSolver(model);
                solver.OrderDofs(false); //v2.1. TODO: Is this needed in this case?
                foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values) linearSystem.Reset();
                //solver.ResetSubdomainForcesVector();
            }


            var elementProvider = new ElementStructuralStiffnessProvider_v2();                      
            #region INTEGRATION constitutive Matrix            
            var integrationSimultaneous = new SubdomainCalculationsAndAssembly();
            (Dictionary<int, double[][]> KfpDqSubdomains, Dictionary<int, double[][]> KppDqVectorsSubdomains) =
                integrationSimultaneous.UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje_v2(model, elementProvider, scaleTransitions, boundaryNodes, boundaryElements, solver);

            Dictionary<int, double[][]> f2_vectorsSubdomains = SubdomainCalculationsMultiple_v2.CalculateKffinverseKfpDqSubdomains_v2(KfpDqSubdomains, model, elementProvider, scaleTransitions, boundaryNodes, solver);

            Dictionary<int, double[][]> f3_vectorsSubdomains = SubdomainCalculationsMultiple_v2.CalculateKpfKffinverseKfpDqSubdomains_v2(f2_vectorsSubdomains, model, elementProvider, scaleTransitions, boundaryNodes);

            double[][] f3_vectors = SubdomainCalculationsMultiple_v2.CombineMultipleSubdomainsIntegrationVectorsIntoTotal_v2(f3_vectorsSubdomains, scaleTransitions);
            double[][] KppDqVectors = SubdomainCalculationsMultiple_v2.CombineMultipleSubdomainsIntegrationVectorsIntoTotal_v2(KppDqVectorsSubdomains, scaleTransitions);

            double[][] f4_vectors = SubdomainCalculations_v2.SubtractConsecutiveVectors_v2(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations_v2.CalculateDqCondDq_v2(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] d2W_dfdf = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    d2W_dfdf[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            
            #endregion

            #region constitutive tensors transformation methods
            double[,] d2W_dFtrdFtr = Reorder_d2Wdfdf_to_d2W_dFtrdFtr(d2W_dfdf);

            double[,] SPK_mat = new double[3, 3] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }; double[,] DefGradMat=new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[,] Cinpk = Transform_d2WdFtrdFtr_to_Cijrs(d2W_dFtrdFtr, SPK_mat, DefGradMat); // to onomazoume Cinpk epeidh einai to 9x9 kai to diakrinoume etsi apo to Cijrs 6x6

            double[,] Cijrs = CombineCinpkTensorTermsIntoMatrix(Cinpk);
           
            constitutiveMatrix = Matrix.CreateFromArray(Cijrs);
            #endregion

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged();

            if (EstimateOnlyLinearResponse)
            {
                model = null;
                boundaryElements = null;
                boundaryNodes = null;
                rveBuilder = null;
                uInitialFreeDOFDisplacementsPerSubdomain = null;
                initialConvergedBoundaryDisplacements = null;
                Cijrs_prev = null;
            }
        }

        #region Print methods for debug
        //private void PrintMethodsForDebug(double[][] KfpDq, double[][] f2_vectors, double[][] f3_vectors, double[][] KppDqVectors, double[][] f4_vectors, double[,] DqCondDq, double[,] d2W_dfdf, double[,] Cijrs)
        //{
        //    string string0 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integration\d2\";

        //    string string1 = String.Concat(string0, @"KfpDq_{0}.txt");

        //    for (int i2 = 0; i2 < KfpDq.GetLength(0); i2++)
        //    {
        //        string path = string.Format(string1, (i2 + 1).ToString());
        //        Vector data = new Vector(KfpDq[i2]);
        //        data.WriteToFile(path);
        //    }

        //    string string2 = String.Concat(string0, @"KffInvKfpDq_{0}.txt");



        //    for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
        //    {
        //        string path = string.Format(string2, (i2 + 1).ToString());
        //        Vector data = new Vector(f2_vectors[i2]);
        //        data.WriteToFile(path);
        //    }

        //    string string3 = String.Concat(string0, @"f3_vectors_{0}.txt");
        //    string string4 = String.Concat(string0, @"KppDqVectors_{0}.txt");
        //    string string5 = String.Concat(string0, @"f4_vectors_{0}.txt");

        //    for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
        //    {
        //        string path = string.Format(string3, (i2 + 1).ToString());
        //        Vector data = new Vector(f3_vectors[i2]);
        //        data.WriteToFile(path);

        //        path = string.Format(string4, (i2 + 1).ToString());
        //        data = new Vector(KppDqVectors[i2]);
        //        data.WriteToFile(path);

        //        path = string.Format(string5, (i2 + 1).ToString());
        //        data = new Vector(f4_vectors[i2]);
        //        data.WriteToFile(path);

        //    }

        //    PrintUtilities.WriteToFile(DqCondDq, String.Concat(string0, @"DqCondDq.txt"));
        //    PrintUtilities.WriteToFile(d2W_dfdf,  String.Concat(string0, @"d2W_dfdf.txt"));
        //    PrintUtilities.WriteToFile(Cijrs, String.Concat(string0, @"Cijrs.txt"));
        //}
        #endregion


    }
    //Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObj_v2
    //Origin  aplo copy apo nl_elements_test.
    //modifications apo UseBase egine UseBaseSimuRand me odhgo Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimuRand se sxesh me to Microstru...Transformation.cs




}
