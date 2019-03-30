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
    /// Primary multiscale analysis class that connects all nesessary structures for a FE2 simulation of shell structures 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class MicrostructureShell2D : StructuralProblemsMicrostructureBase_v2, IShellMaterial_v2 //A.1
    {
        public double[] NormalVectorV3 { get; set; }
        public double[] TangentVectorV1 { get; set; }
        public double[] TangentVectorV2 { get; set; }

        private Model_v2 model { get; set; }
        private Dictionary<int, Node_v2> boundaryNodes { get; set; }
        Dictionary<int, Dictionary<int, Element_v2>> boundaryElements;
        private IdegenerateRVEbuilder_v2 rveBuilder;
        private bool EstimateOnlyLinearResponse;
        private double volume;
        Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain;
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        private IScaleTransitions_v2 scaleTransitions = new SmallStrain3Dto2DplaneStressScaleTransition(); //TODO: mporoume na to dinoume ston constructor
        Random rnd1 = new Random();
        private readonly Func<Model_v2, ISolver_v2> createSolver;

        Matrix constitutiveMatrix;
        private double[] trueStressVec; // TODO: rename stresses 
        Matrix transformationMatrix; // gia to shell
        private bool modified; // opws sto MohrCoulomb gia to modified

        private double[,] Cijrs_prev;
        private bool matrices_not_initialized = true;
        private double tol;
        public void InitializeMatrices()
        {
            Cijrs_prev = new double[3,3];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
            constitutiveMatrix = Matrix.CreateZero(3,3);
            this.CalculateTransformationMatrix(Vector.CreateFromArray(TangentVectorV1),Vector.CreateFromArray(TangentVectorV2));
        }

        
        //Random properties 
        private int database_size;

        public MicrostructureShell2D(IdegenerateRVEbuilder_v2 rveBuilder, Func<Model_v2, ISolver_v2> createSolver, 
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
                var RigidBodyNodeConstraints = rveBuilder.GetModelRigidBodyNodeConstraints(model);
                scaleTransitions.ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(model, boundaryNode, RigidBodyNodeConstraints);
            }
        }

        private void InitializeFreeAndPrescribedDofsInitialDisplacementVectors()
        {
            uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, IVector>();
            foreach(Subdomain_v2 subdomain in model.SubdomainsDictionary.Values)
            {
                uInitialFreeDOFDisplacementsPerSubdomain.Add(subdomain.ID, Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs));// prosoxh sto Id twn subdomain
            }            
            double[] smallStrainVec = new double[3] ;
            initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, initialConvergedBoundaryDisplacements);
            }            
        }

        public IShellMaterial_v2 Clone()
        {
            Random rnd1 = new Random();
            int new_rve_id = rnd1.Next(1, database_size+1);
            return new MicrostructureShell2D((IdegenerateRVEbuilder_v2)rveBuilder.Clone(new_rve_id), createSolver, EstimateOnlyLinearResponse,database_size);
        }

        object ICloneable.Clone() => this.Clone();

        public Dictionary<int, Node_v2> BoundaryNodesDictionary
        {
            get { return boundaryNodes; }
        }
        public IList<Node_v2> BoundaryNodes
        {
            get { return boundaryNodes.Values.ToList<Node_v2>(); }
        }

        public void UpdateMaterial(double[] smallStrainVec)
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

            double[] rveCoordinatesSmallStrainVec = TransformStrains(smallStrainVec);
            smallStrainVec = rveCoordinatesSmallStrainVec;

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    Cijrs_prev[i1, j1] = this.constitutiveMatrix[i1, j1];}
            }

            #region Rve prescribed Dofs total DIsplacement Dictionary Creation (nessesary for NRNLAnalyzer)
            Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node_v2 boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, totalPrescribedBoundaryDisplacements);
            }
            #endregion


            //var linearSystems = CreateNecessaryLinearSystems(model);          
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

            trueStressVec = new double [DqFpp.Length];
            for (int i1 = 0; i1 < DqFpp.Length; i1++)
            { trueStressVec[i1]=(1 / volume) * DqFpp[i1]; }

            #endregion

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
            // transformation gia to shell 
            (var transformedTrueStressVec, var transformedConstitutiveMat) =StressesAndConstitutiveMatrixTransformation(trueStressVec, constitutiveMat);
            #endregion

            this.constitutiveMatrix = Matrix.CreateFromArray(transformedConstitutiveMat);
            trueStressVec = transformedTrueStressVec;

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged(); 
        }

        
        private void CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = Matrix.CreateZero(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            Matrix inverse = auxMatrix1.Invert(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            Ei[0] = surfaceBasisVector1.CopyToArray();
            double G1_norm = surfaceBasisVector1.Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = (Vector.CreateFromArray(Ei[1])).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = Vector.CreateFromArray(Ei[i1]).DotProduct(Vector.CreateFromArray(G_i[i2]));
                }
            }

             transformationMatrix = Matrix.CreateFromArray(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });
        }

        private double[] TransformStrains(double[] smallStrainVec)
        {
            double[] rveCoordinatesSmallStrainVec = (transformationMatrix *(Vector.CreateFromArray( smallStrainVec))).CopyToArray();
            return rveCoordinatesSmallStrainVec;
        }


        private (double[] , double[,] ) StressesAndConstitutiveMatrixTransformation( double[] trueStressVec, double[,] constitutiveMat)
        {            
            var transformedTrueStressVec = (transformationMatrix.Transpose() * (Vector.CreateFromArray(trueStressVec))).CopyToArray();
            // TODO: CHECK matrix for corrupted data after multiplication
            var transformedConstitutiveMat = (transformationMatrix.Transpose() * (Matrix.CreateFromArray(constitutiveMat)) * transformationMatrix).CopyToArray2D();

            return (transformedTrueStressVec, transformedConstitutiveMat);
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
                return constitutiveMatrix; // TODO: apla kratame to constitutive matrix san array[,] (alla matrix mporei na xrhimopoithei gia tis peristrofes)
            }
        }

        public double[] Stresses // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return  trueStressVec; }
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

        public double YoungModulus => throw new NotImplementedException(); // TODO: remove these from the interface.

        public double PoissonRatio => throw new NotImplementedException();


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

        public void CalculateOriginalConstitutiveMatrixWithoutNLAnalysis()
        {
            ISolver_v2 solver;
            if (matrices_not_initialized)
            {
                this.InitializeMatrices();
                this.InitializeData();
                solver = createSolver(model);
                solver.OrderDofs(false); //model.GlobalDofOrdering = solver.DofOrderer.OrderDofs(model); //TODO find out if new structures cause any problems
                foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
                {
                    linearSystem.Reset();
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

            trueStressVec = new double[3];

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
            // transformation gia to shell 
            (var transformedTrueStressVec, var transformedConstitutiveMat) = StressesAndConstitutiveMatrixTransformation(trueStressVec, constitutiveMat);
            this.constitutiveMatrix = Matrix.CreateFromArray(transformedConstitutiveMat);
            trueStressVec = transformedTrueStressVec;
            #endregion

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged();

            if (EstimateOnlyLinearResponse)
            {
                model = null;
                boundaryElements = null;
                boundaryNodes = null;
                NormalVectorV3 = null;
                TangentVectorV1 = null;
                TangentVectorV2 = null;
                rveBuilder = null;
                uInitialFreeDOFDisplacementsPerSubdomain = null;
                initialConvergedBoundaryDisplacements = null;
                trueStressVec = null;
                transformationMatrix = null;
                Cijrs_prev = null;
            }

        }

        //todo delete when unnesessary
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
    //Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimuRandObj_v2
    // Origin: Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimuRand
    // modifications: --> updated se v2



}
