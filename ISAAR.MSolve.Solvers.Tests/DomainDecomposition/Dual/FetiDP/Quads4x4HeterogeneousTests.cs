using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4HeterogeneousTests
    {
        [Fact]
        public static void TestScalingMatrices()
        {
            #region Replace the next with hardcoded matrices and mocking objects
            // Run the analysis so that all objects are created
            // Setup the model
            double stiffnessRatio = 1E-2; // Do not change this! The expected solution is taken for this value
            Model model = CreateModel(stiffnessRatio);
            Dictionary<int, HashSet<INode>> cornerNodes = Quads4x4MappingMatricesTests.DefineCornerNodes(model);

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
            var fetiSolverBuilder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            var preconditionerFactory = new DirichletPreconditioner.Factory();
            fetiSolverBuilder.PreconditionerFactory = preconditionerFactory;
            FetiDPSolver fetiSolver = fetiSolverBuilder.BuildSolver(model);

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();
            #endregion

            // Access private fields of FetiDPSolver and DirichletPreconditioner.Factory using reflection
            FieldInfo fi = typeof(FetiDPSolver).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumerator)fi.GetValue(fetiSolver);
            fi = typeof(FetiDPSolver).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparator)fi.GetValue(fetiSolver);
            fi = typeof(FetiDPSolver).GetField("stiffnessDistribution", BindingFlags.NonPublic | BindingFlags.Instance);
            var stiffnessDistribution = (IStiffnessDistribution)fi.GetValue(fetiSolver);
            MethodInfo method = preconditionerFactory.GetType().GetMethod("CalcBoundaryPreconditioningBooleanMatrices",
                BindingFlags.NonPublic | BindingFlags.Instance);
            var Bpbr = (Dictionary<int, IMappingMatrix>)method.Invoke(preconditionerFactory, 
                new object[] { stiffnessDistribution, dofSeparator, lagrangeEnumerator });

            // Compare the mapping matrices against the expected ones
            var expectedBpbr = new Dictionary<int, Matrix>();
            expectedBpbr[0] = Matrix.CreateFromArray(new double[,]
            {
                { 0.5, 0, 0, 0 },
                { 0, 0.5, 0, 0 },
                { 0, 0, 0.5, 0 },
                { 0, 0, 0, 0.5 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 }
            });
            expectedBpbr[1] = Matrix.CreateFromArray(new double[,]
            {
                { -0.5, 0, 0, 0 },
                { 0, -0.5, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0.00990099009900990, 0 },
                { 0, 0, 0, 0.00990099009900990 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 }
            });
            expectedBpbr[2] = Matrix.CreateFromArray(new double[,]
            {
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { -0.5, 0, 0, 0 },
                { 0, -0.5, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0.00990099009900990, 0 },
                { 0, 0, 0, 0.00990099009900990 }
            });
            expectedBpbr[3] = Matrix.CreateFromArray(new double[,]
            {
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { -0.990099009900990, 0, 0, 0 },
                { 0, -0.990099009900990, 0, 0 },
                { 0, 0, -0.990099009900990, 0 },
                { 0, 0, 0, -0.990099009900990 }
            });

            double tol = 1E-13;
            foreach (int s in expectedBpbr.Keys)
            {
                Matrix explicitBpr = Bpbr[s].MultiplyRight(Matrix.CreateIdentity(Bpbr[s].NumColumns));
                Assert.True(expectedBpbr[s].Equals(explicitBpr, tol));
            }
        }

        [Fact]
        public static void TestSolver()
        {
            // Setup the model
            double stiffnessRatio = 1E-2; // Do not change this! The expected solution is taken for this value
            Model model = CreateModel(stiffnessRatio);
            Dictionary<int, HashSet<INode>> cornerNodes = Quads4x4MappingMatricesTests.DefineCornerNodes(model);

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
            var fetiSolverBuilder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();
            FetiDPSolver fetiSolver = fetiSolverBuilder.BuildSolver(model);

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);

            // Check against expected solution
            double tol = 1E-7;
            var globalUExpected = Vector.CreateFromArray(new double[]
            {
                17.623494584618864, 12.564560593215612, 31.832863897135404, 34.496634608059082, 40.255481382985629,
                66.49190654178912, 42.572002358887204, 99.798764204232072, 4.267568672307144, 9.00506902466324,
                9.100928263505315, 31.107370029452451, 12.1615036308774, 66.065492717632239, 11.510673148931499,
                102.06649895017948, -3.0529124682202156, 9.24107474483673, -7.8531777412741217, 26.728892403726846,
                -16.890006178831449, 70.602493468916791, -21.80233265288679, 109.39882637058051, -4.7311061272016808,
                10.030926199331375, -5.6722429958962142, 18.837815470700932, 146.94209278892487, 392.04674590737193,
                -35.619167413693908, 1407.200332011206, -9.9609496807814057, 10.46574373452243, -17.603838651152756,
                20.760800663270086, -843.13592713307355, 371.10700308359418, -1666.2547486301742, 3714.1637893447919
            });
            Assert.True(globalUExpected.Equals(globalU, tol));
        }

        private static Model CreateModel(double stiffnessRatio)
        {
            //                                    Λ P
            //                                    | 
            //                                     
            // |> 20 ---- 21 ---- 22 ---- 23 ---- 24
            //    |  (12) |  (13) |  (14) |  (15) |
            //    |  E0   |  E0   |  E1   |  E1   |
            // |> 15 ---- 16 ---- 17 ---- 18 ---- 19
            //    |  (8)  |  (9)  |  (10) |  (11) |
            //    |  E0   |  E0   |  E1   |  E1   |
            // |> 10 ---- 11 ---- 12 ---- 13 ---- 14
            //    |  (4)  |  (5)  |  (6)  |  (7)  |
            //    |  E0   |  E0   |  E0   |  E0   |
            // |> 5 ----- 6 ----- 7 ----- 8 ----- 9
            //    |  (0)  |  (1)  |  (2)  |  (3)  |
            //    |  E0   |  E0   |  E0   |  E0   |
            // |> 0 ----- 1 ----- 2 ----- 3 ----- 4


            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 4.0;
            builder.DomainLengthY = 4.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 4;
            builder.NumTotalElementsY = 4;
            //builder.YoungModulus = 1.0;
            double E0 = 1.0;
            builder.YoungModuliOfSubdomains = new double[,]
            {
                { E0, E0 }, { E0, stiffnessRatio * E0}
            };
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.UpperRightCorner, StructuralDof.TranslationY, 10.0);

            return builder.BuildModel();
        }
    }
}
