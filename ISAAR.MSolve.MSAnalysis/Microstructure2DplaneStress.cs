using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Primary multiscale analysis class that connects all nesessary structures for a FE2 simulation (with an 3D to 2D degenerate Rve)
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class Microstructure2DplaneStress : StructuralProblemsMicrostructureBase, IContinuumMaterial2D
    {
        private Model model { get; set; }
        //private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private Dictionary<int, Node> boundaryNodes { get; set; }
        Dictionary<int, Dictionary<int, Element>> boundaryElements;
        private IdegenerateRVEbuilder rveBuilder;
        private bool EstimateOnlyLinearResponse;
        //private NewtonRaphsonNonLinearAnalyzer microAnalyzer;
        private double volume;
        public Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain { get; private set; }
        Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements;
        private IScaleTransitions scaleTransitions = new SmallStrain3Dto2DplaneStressScaleTransition();
        Random rnd1 = new Random();
        private readonly Func<Model, ISolver> createSolver;

        // aparaithta gia to implementation tou IFiniteElementMaterial3D
        Matrix constitutiveMatrix;
        private double[] trueStressVec = new double[6];
        private bool modified; // opws sto MohrCoulomb gia to modified

        private double[,] Cijrs_prev;
        private bool matrices_not_initialized = true;
        private double tol;
        public void InitializeMatrices()
        {
            Cijrs_prev = new double[3, 3];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
            constitutiveMatrix = Matrix.CreateZero(3,3);
        }


        //double[] Stresses { get; }
        //IMatrix2D ConstitutiveMatrix { get; } TODOGerasimos

        //Random properties 
        private int database_size;

        public Microstructure2DplaneStress(IdegenerateRVEbuilder rveBuilder, Func<Model, ISolver> createSolver, bool EstimateOnlyLinearResponse, int database_size)
        {
            this.rveBuilder = rveBuilder;
            this.createSolver = createSolver;
            this.EstimateOnlyLinearResponse = EstimateOnlyLinearResponse;
            this.database_size = database_size;            
        }

        private void InitializeData()
        {
            Tuple<Model, Dictionary<int, Node>, double> modelAndBoundaryNodes = this.rveBuilder.GetModelAndBoundaryNodes();
            this.model = modelAndBoundaryNodes.Item1;
            this.boundaryNodes = modelAndBoundaryNodes.Item2;
            this.boundaryElements = GetSubdomainsBoundaryFiniteElementsDictionaries(model, boundaryNodes);
            this.volume = modelAndBoundaryNodes.Item3;
            DefineAppropriateConstraintsForBoundaryNodes();
            this.model.ConnectDataStructures();            
        }


        private void DefineAppropriateConstraintsForBoundaryNodes()
        {
            var RigidBodyNodeConstraints = rveBuilder.GetModelRigidBodyNodeConstraints(model);
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(model, boundaryNode, RigidBodyNodeConstraints);
            }
        }

        private void InitializeFreeAndPrescribedDofsInitialDisplacementVectors()
        {
            uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, IVector>();
            foreach(Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                uInitialFreeDOFDisplacementsPerSubdomain.Add(subdomain.ID, Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs));// prosoxh sto Id twn subdomain
            }
            double[] smallStrainVec = new double[3];
            initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<IDofType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, initialConvergedBoundaryDisplacements);
            }            
        }

        public object Clone()
        {
            int new_rve_id = rnd1.Next(1, database_size + 1);
            return new Microstructure2DplaneStress((IdegenerateRVEbuilder)rveBuilder.Clone(new_rve_id), createSolver, EstimateOnlyLinearResponse, database_size);
        }

        public Dictionary<int, Node> BoundaryNodesDictionary
        {
            get { return boundaryNodes; }
        }
        public IList<Node> BoundaryNodes
        {
            get { return boundaryNodes.Values.ToList<Node>(); }
        }

        public void UpdateMaterial(double[] smallStrainVec)
        {
            ISolver solver;
            if (matrices_not_initialized)
            {
                this.InitializeMatrices();
                this.InitializeData();
                solver = createSolver(model);
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in solver.LinearSystems.Values)
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
                //solver.ResetSubdomainForcesVector();
                foreach (ILinearSystem linearSystem in solver.LinearSystems.Values)
                {
                    linearSystem.Reset();
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces; //TODO MS 
                }
            }

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {Cijrs_prev[i1, j1] = constitutiveMatrix[i1, j1];}
            }

            #region Rve prescribed Dofs total DIsplacement Dictionary Creation (nessesary for NRNLAnalyzer)
            // epivolh metakinhsewn ston analyzer pou molis dhmiourghsame --> tha ginetai sto dictionary
            Dictionary<int, Dictionary<IDofType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<IDofType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, totalPrescribedBoundaryDisplacements);
            }
            #endregion
                     
            //var linearSystems = CreateNecessaryLinearSystems(model);    // OPOU pairnei rhs apo subdomainForces       
            //var solver = GetAppropriateSolver(linearSystems);

            
            #region Creation of nessesary analyzers for NRNLAnalyzer and Creation of Microstructure analyzer (NRNLdevelop temporarilly) and solution ;
            int increments = 1; int MaxIterations = 100; int IterationsForMatrixRebuild = 1;
            (MicrostructureBvpNRNLAnalyzer microAnalyzer, ProblemStructural provider, ElementStructuralStiffnessProvider elementProvider) = 
                AnalyzeMicrostructure(model, solver, increments, MaxIterations, IterationsForMatrixRebuild,
                totalPrescribedBoundaryDisplacements, initialConvergedBoundaryDisplacements, boundaryNodes, uInitialFreeDOFDisplacementsPerSubdomain);
            
            #endregion

            #region update of free converged displacements vectors
            uInitialFreeDOFDisplacementsPerSubdomain = microAnalyzer.GetConvergedSolutionVectorsOfFreeDofs();// ousiastika to u pou twra taftizetai me to uPlusuu
            #endregion


            #region INTEGRATION stresses 
            Dictionary<int, IVector> du = microAnalyzer.GetConvergedIncrementalSolutionVectorsOfFreeDofs();
            Dictionary<int, double[]> FppReactionVectorSubdomains = SubdomainCalculationsMultiple.CalculateFppReactionsVectorSubdomains(model, elementProvider, scaleTransitions, boundaryNodes,
                uInitialFreeDOFDisplacementsPerSubdomain, du, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, increments, increments);
            double[] FppReactionVector= SubdomainCalculationsMultiple.CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal(FppReactionVectorSubdomains);



            double[] DqFpp = SubdomainCalculations.CalculateDqFpp(FppReactionVector, scaleTransitions, boundaryNodes);

            trueStressVec = new double [DqFpp.Length];
            for (int i1 = 0; i1 < DqFpp.Length; i1++)
            { trueStressVec[i1]=(1 / volume) * DqFpp[i1]; }

           
            //SPK_vec = new double[6] { SPK_mat[0,0], SPK_mat[1,1], SPK_mat[2,2], SPK_mat[0,1], SPK_mat[1,2], SPK_mat[0,2] };
            //TODOna elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            #endregion

            #region INTEGRATION constitutive Matrix
            var integrationSimultaneous = new SubdomainCalculationsAndAssembly();
            (Dictionary<int, double[][]> KfpDqSubdomains, Dictionary<int, double[][]> KppDqVectorsSubdomains) = 
                integrationSimultaneous.UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje(model, elementProvider, scaleTransitions, boundaryNodes, boundaryElements, solver);


            Dictionary<int, double[][]> f2_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKffinverseKfpDqSubdomains(KfpDqSubdomains, model, elementProvider, scaleTransitions, boundaryNodes, solver);

            Dictionary<int, double[][]> f3_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKpfKffinverseKfpDqSubdomains(f2_vectorsSubdomains, model, elementProvider, scaleTransitions, boundaryNodes);

            double[][] f3_vectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(f3_vectorsSubdomains,scaleTransitions);
            double[][] KppDqVectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(KppDqVectorsSubdomains,scaleTransitions);

            double[][] f4_vectors = SubdomainCalculations.SubtractConsecutiveVectors(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations.CalculateDqCondDq(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] constitutiveMat = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    constitutiveMat[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            #endregion

            #region update of prescribed converged displacements vectors;
            initialConvergedBoundaryDisplacements = totalPrescribedBoundaryDisplacements;
            #endregion

            #region constitutive tensors transformation methods
            
            
            #endregion

            this.constitutiveMatrix = Matrix.CreateFromArray(constitutiveMat);

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged(); 
        }        

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
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
            get { return trueStressVec; }
        }

        public void SaveState()
        {
            var subdomainUpdaters = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions>(1); 
            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            {
                subdomainUpdaters.Add(subdomain.ID, new NonLinearSubdomainUpdaterWithInitialConditions(subdomain)); //v2.3
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


        

        public void CalculateOriginalConstitutiveMatrixWithoutNLAnalysis()
        {
            ISolver solver;
            if (matrices_not_initialized)
            {
                this.InitializeMatrices();
                this.InitializeData();
                solver = createSolver(model);
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in solver.LinearSystems.Values)
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
                foreach (ILinearSystem linearSystem in solver.LinearSystems.Values) linearSystem.Reset();
                //solver.ResetSubdomainForcesVector();
            }


            var elementProvider = new ElementStructuralStiffnessProvider();                      
            #region INTEGRATION constitutive Matrix            
            var integrationSimultaneous = new SubdomainCalculationsAndAssembly();
            (Dictionary<int, double[][]> KfpDqSubdomains, Dictionary<int, double[][]> KppDqVectorsSubdomains) =
                integrationSimultaneous.UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje(model, elementProvider, scaleTransitions, boundaryNodes, boundaryElements, solver);

            Dictionary<int, double[][]> f2_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKffinverseKfpDqSubdomains(KfpDqSubdomains, model, elementProvider, scaleTransitions, boundaryNodes, solver);

            Dictionary<int, double[][]> f3_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKpfKffinverseKfpDqSubdomains(f2_vectorsSubdomains, model, elementProvider, scaleTransitions, boundaryNodes);

            double[][] f3_vectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(f3_vectorsSubdomains, scaleTransitions);
            double[][] KppDqVectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(KppDqVectorsSubdomains, scaleTransitions);

            double[][] f4_vectors = SubdomainCalculations.SubtractConsecutiveVectors(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations.CalculateDqCondDq(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] constitutiveMat = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    constitutiveMat[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            
            #endregion

            #region constitutive tensors transformation methods                       
            constitutiveMatrix = Matrix.CreateFromArray(constitutiveMat);
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

        #region transformation methods
        //TODO: implement and use methods for shell transformation
        private double[] transformTrueStressVec(double[] trueStressVec, double[] tangent1, double[] tangent2, double[] normal)
        {
            throw new NotImplementedException();
        }

        private double[,] TransformConstitutiveMatrix(double[,] constitutiveMat, double[] tangent1, double[] tangent2, double[] normal)
        {
            throw new NotImplementedException();
        }

        #endregion

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
    //Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObjSmallStrains2DplaneStress
    //Origin  Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObj
    //modifications apo defgrad egine smallstrains2d me odhgo Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains2D se sxesh me to Microstructure3DevelopMultipleSubdomainsUseBase.cs




}
