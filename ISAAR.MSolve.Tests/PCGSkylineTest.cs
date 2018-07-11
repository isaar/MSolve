using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.PCGSkyline;
using ISAAR.MSolve.Solvers.PCG;
using Moq;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Tests
{
    public class PCGSkylineTest
    {
        [Fact]
        public void TestPCGSkylineExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            int totalDOFs = 50;
            int totalNodes = totalDOFs / 2;
            int totalElements = totalNodes - 1;
            int iNodeID = 0;
            int iElementID = 0;
            for (int iNode = 0; iNode < totalNodes; iNode++)
            {
                iNodeID++;
                var n = new Node() { ID = iNodeID, X = iNodeID * 10, Y = 0, Z = 0 };
                model.NodesDictionary.Add(iNodeID, n);
            }
            iNodeID = 0;
            for (int iElement = 0; iElement < totalElements; iElement++)
            {
                iElementID++;
                iNodeID++;
                var e = new Element() { ID = iElementID };
                e.NodesDictionary.Add(iNodeID, model.NodesDictionary[iNodeID]);
                e.NodesDictionary.Add(iNodeID + 1, model.NodesDictionary[iNodeID + 1]);
                var m = new Mock<IFiniteElement>();
                //---------------------------------------------------------------------------------
                m.Setup(x => x.StiffnessMatrix(e)).Returns(new SkylineMatrix2D(new double[,] {
                    { 2, -1, 0, 0 }, { -1, 2, -1, 0 }, { 0, -1, 2, -1 }, { 0, 0, -1, 2 }
                }));
                m.Setup(x => x.GetElementDOFTypes(e)).Returns(new[] { new[] { DOFType.X, DOFType.Y }, new[] { DOFType.X, DOFType.Y } });
                m.SetupGet(x => x.DOFEnumerator).Returns(new GenericDOFEnumerator());
                e.ElementType = m.Object;
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            }

            var linearSystems = new Dictionary<int, ILinearSystem>();
            double[] externalForces = new double[totalDOFs];
            externalForces[totalDOFs - 1] = +10;
            linearSystems.Add(1, new SkylineLinearSystem(1, externalForces));
            model.ConnectDataStructures();
            SolverPCGSimpleSearchVectorCalculator search = new SolverPCGSimpleSearchVectorCalculator();
            SolverPCG<Matrix2D> solver = new SolverPCG<Matrix2D>(linearSystems[1], search);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer staticAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);
            staticAnalyzer.BuildMatrices();
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();
            var v = linearSystems[1].Solution;
            var b = new double[50];
            linearSystems[1].Matrix.Multiply(v, b);
            var r = ((Vector)linearSystems[1].RHS).Data;
            var solutionTolerance = 1e-6;
            var elementsOutsideTolerance = b.Zip(r, (valueB, valueR) => Math.Abs(valueB - valueR)).Count(x => x > solutionTolerance);
            Assert.Equal(elementsOutsideTolerance, 0);
        }
    }
}

