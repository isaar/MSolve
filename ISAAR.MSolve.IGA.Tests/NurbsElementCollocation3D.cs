using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Assemblers.Collocation;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsElementCollocation3D
    {
        [Fact]
        public void CollocationPoint3DMatrix()
        {
            var model = new CollocationModel();
            ModelCreator modelCreator = new ModelCreator(model);
            string filename = "..\\..\\..\\InputFiles\\Collocation 3D.txt";
            IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
            modelReader.CreateCollocationModelFromFile();

            var gmresBuilder =new GmresSolver.Builder();
            ISolver solver = gmresBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.BuildMatrices();

            var k = solver.LinearSystems[0].Matrix;
            Matrix<double> stiffnessMatrixExpected =
                MatlabReader.Read<double>("..\\..\\..\\InputFiles\\Kcol3D.mat", "Ktotal");
            

            for (int i = 0; i < k.NumRows; i++)
            {
                for (int j = 0; j < k.NumColumns; j++)
                {
                    Utilities.AreValuesEqual(stiffnessMatrixExpected[i, j], k[i, j], 10e-9);
                }
            }

        }


        [Fact]
        public void GmresTest()
        {
            IMatrixView matrix = Matrix.CreateFromArray(new double[9,9]
            {
                { 2,  0,  0, -1,  0,  0,  0,  0,  0},
                { 0,  2, -1,  0,  0,  0,  0,  0,  0},
                { 0, -1,  2,  0,  0,  0,  0,  0,  0},
                { -1,  0,  0,  2, -1,  0,  0,  0,  0},
                { 0,  0,  0, -1,  2, -1,  0,  0,  0},
                { 0,  0,  0,  0, -1,  2, -1,  0,  0},
                { 0,  0,  0,  0,  0, -1,  2, -1,  0},
                { 0,  0,  0,  0,  0,  0, -1,  2, -1},
                { 0,  0,  0,  0,  0,  0,  0, -1,  2}
            });
            IVectorView rhs = Vector.CreateWithValue(9, 1.0);
            var solution = Vector.CreateZero(9);
            var gmresAlgorithmBuilder= new GmresAlgorithm.Builder()
            {
                MaximumIterations = 20,
                InnerIterationsProvider = new FixedMaxIterationsProvider(8),
                AbsoluteTolerance = 1e-8,
                RelativeTolerance = 1e-8
            };
            var gmres = gmresAlgorithmBuilder.Build();
            gmres.Solve(matrix, rhs, solution,true, () => Vector.CreateZero(9));

            var expectedSolution = Vector.CreateFromArray(new double[] {3.5, 1.0, 1.0, 6.0, 7.5, 8.0, 7.5, 6.0, 3.5});

            for (int i = 0; i < expectedSolution.Length; i++)
            {
                Utilities.AreValuesEqual(expectedSolution[i], solution[i], 10e-9);
            }

        }
    }
}