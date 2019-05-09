using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Readers;
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
            Model model = new Model();
            ModelCreator modelCreator = new ModelCreator(model);
            string filename = "..\\..\\..\\InputFiles\\Collocation 3D.txt";
            IsogeometricReader modelReader = new IsogeometricReader(modelCreator, filename);
            modelReader.CreateCollocationModelFromFile();

            ISolver solver = new GmresSolver(model,
                new AsymmetricDofOrderer(new RowDofOrderingStrategy()),
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()),
                new CsrRectangularAssembler(), "CsrRectangularAssembler");

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
    }
}