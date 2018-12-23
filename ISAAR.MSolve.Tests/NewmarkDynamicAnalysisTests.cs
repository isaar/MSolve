using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Moq;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using Xunit;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.Tests
{
    public static class NewmarkDynamicAnalysisTests
    {
        [Fact]
        public static void TestBatheImplicitAnalysisExample()
        {
            Numerical.LinearAlgebra.VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            var n = new Node() { ID = 0 };
            var e = new Element() { ID = 0 };
            e.NodesDictionary.Add(0, n);
            var m = new Mock<IFiniteElement>();
            m.Setup(x => x.StiffnessMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 6, -2, 4 }));
            m.Setup(x => x.MassMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 }));
            m.Setup(x => x.DampingMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 0, 0, 0 }));
            m.Setup(x => x.GetElementDOFTypes(e)).Returns(new[] { new[] { DOFType.X, DOFType.Y } });
            m.SetupGet(x => x.DOFEnumerator).Returns(new GenericDOFEnumerator());
            e.ElementType = m.Object;
            model.NodesDictionary.Add(0, n);
            model.ElementsDictionary.Add(0, e);
            model.SubdomainsDictionary[0].ElementsDictionary.Add(0, e);
            model.Loads.Add(new Load() { Amount = 10, Node = n, DOF = DOFType.Y });
            var lX = new Mock<IMassAccelerationHistoryLoad>();
            lX.SetupGet(x => x.DOF).Returns(DOFType.X);
            lX.SetupGet(x => x[It.IsAny<int>()]).Returns(0);
            var lY = new Mock<IMassAccelerationHistoryLoad>();
            lY.SetupGet(x => x.DOF).Returns(DOFType.Y);
            lY.SetupGet(x => x[0]).Returns(10);
            lY.SetupGet(x => x[It.IsInRange(1, 100, Range.Inclusive)]).Returns(0);
            model.MassAccelerationHistoryLoads.Add(lX.Object);
            model.MassAccelerationHistoryLoads.Add(lY.Object);
            model.ConnectDataStructures();
            m.Setup(x => x.CalculateAccelerationForces(It.IsAny<Element>(), It.IsAny<IList<MassAccelerationLoad>>()))
                .Returns<Element, IList<MassAccelerationLoad>>((element, loads) =>
                {
                    Numerical.LinearAlgebra.Vector accelerations = new Numerical.LinearAlgebra.Vector(2);
                    accelerations[0] = loads[0].Amount;
                    accelerations[1] = loads[1].Amount;

                    var massMatrix = new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 });
                    double[] forces = new double[2];
                    massMatrix.Multiply(accelerations, forces);
                    return forces;
                }
            );

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems.Add(0, new SkylineLinearSystem(0, new[] { 0.0, 10.0 }));
            linearSystems[0].Matrix = new Numerical.LinearAlgebra.SkylineMatrix2D(new double[,] { { 6, -2 }, { -2, 4 } });

            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer dynamicAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            dynamicAnalyzer.BuildMatrices();
            dynamicAnalyzer.Initialize();
            dynamicAnalyzer.Solve();
            Assert.Equal(2.2840249264795207, linearSystems[0].Solution[0], 8);
            Assert.Equal(2.4351921891904156, linearSystems[0].Solution[1], 8);
        }

        [Fact]
        private static void TestBatheImplicitAnalysisExample_v2()
        {
            var model = new Model_v2();
            int subdomainID = 0;
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            var n = new Node() { ID = 0 };
            var e = new Element() { ID = 0 };
            e.NodesDictionary.Add(0, n);
            var m = new Mock<IFiniteElement>();
            m.Setup(x => x.StiffnessMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 6, -2, 4 }));
            m.Setup(x => x.MassMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 }));
            m.Setup(x => x.DampingMatrix(e)).Returns(new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 0, 0, 0 }));
            m.Setup(x => x.GetElementDOFTypes(e)).Returns(new[] { new[] { DOFType.X, DOFType.Y } });
            m.SetupGet(x => x.DOFEnumerator).Returns(new GenericDOFEnumerator());
            e.ElementType = m.Object;
            model.NodesDictionary.Add(0, n);
            model.ElementsDictionary.Add(0, e);
            model.SubdomainsDictionary[subdomainID].Elements.Add(e);
            model.Loads.Add(new Load() { Amount = 10, Node = n, DOF = DOFType.Y });
            var lX = new Mock<IMassAccelerationHistoryLoad>();
            lX.SetupGet(x => x.DOF).Returns(DOFType.X);
            lX.SetupGet(x => x[It.IsAny<int>()]).Returns(0);
            var lY = new Mock<IMassAccelerationHistoryLoad>();
            lY.SetupGet(x => x.DOF).Returns(DOFType.Y);
            lY.SetupGet(x => x[0]).Returns(10);
            lY.SetupGet(x => x[It.IsInRange(1, 100, Range.Inclusive)]).Returns(0);
            model.MassAccelerationHistoryLoads.Add(lX.Object);
            model.MassAccelerationHistoryLoads.Add(lY.Object);
            m.Setup(x => x.CalculateAccelerationForces(It.IsAny<Element>(), It.IsAny<IList<MassAccelerationLoad>>()))
                .Returns<Element, IList<MassAccelerationLoad>>((element, loads) =>
                {
                    var accelerations = new Numerical.LinearAlgebra.Vector(2);
                    accelerations[0] = loads[0].Amount;
                    accelerations[1] = loads[1].Amount;

                    var massMatrix = new Numerical.LinearAlgebra.SymmetricMatrix2D(new double[] { 2, 0, 1 });
                    double[] forces = new double[2];
                    massMatrix.Multiply(accelerations, forces);
                    return forces;
                }
            );

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            SkylineSolver solver = solverBuilder.BuildSolver(model);

            solver.LinearSystems[subdomainID].SetMatrix(
                SkylineMatrix.CreateFromArrays(2, new double[] { 6, 4, -2 }, new int[] { 0, 1, 3 }, true)); // K = [6 -2; -2 4]
            solver.LinearSystems[subdomainID].RhsVector = Vector.CreateFromArray(new double[] { 0, 10 });

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer_v2(solver);

            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer_v2.Builder(model, solver, provider, childAnalyzer, 0.28, 3.36);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            //parentAnalyzerBuilder.SetNewmarkParameters(0.25, 0.5); // Not necessary. This is the default
            NewmarkDynamicAnalyzer_v2 parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { 0, 1 });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(2.2840249264795207, log.DOFValues[0], 8);
            Assert.Equal(2.4351921891904156, log.DOFValues[1], 8);
        }
    }
}

