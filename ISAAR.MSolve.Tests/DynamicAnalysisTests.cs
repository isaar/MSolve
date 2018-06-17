using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
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

namespace ISAAR.MSolve.Tests
{
    public class DynamicAnalysisTests
    {
        [Fact]
        public void TestBatheImplicitAnalysisExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            var n = new Node() { ID = 0 };
            var e = new Element() { ID = 0 };
            e.NodesDictionary.Add(0, n);
            var m = new Mock<IFiniteElement>();
            m.Setup(x => x.StiffnessMatrix(e)).Returns(new SymmetricMatrix2D(new double[] { 6, -2, 4 }));
            m.Setup(x => x.MassMatrix(e)).Returns(new SymmetricMatrix2D(new double[] { 2, 0, 1 }));
            m.Setup(x => x.DampingMatrix(e)).Returns(new SymmetricMatrix2D(new double[] { 0, 0, 0 }));
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
                    Vector accelerations = new Vector(2);
                    accelerations[0] = loads[0].Amount;
                    accelerations[1] = loads[1].Amount;

                    var massMatrix = new SymmetricMatrix2D(new double[] { 2, 0, 1 });
                    double[] forces = new double[2];
                    massMatrix.Multiply(accelerations, forces);
                    return forces;
                }
            );


            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems.Add(0, new SkylineLinearSystem(0, new[] { 0, 10d }));
            linearSystems[0].Matrix = new SkylineMatrix2D(new double[,] { { 6, -2 }, { -2, 4 } });
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer dynamicAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);
            //StaticAnalyzer staticAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            dynamicAnalyzer.BuildMatrices();
            dynamicAnalyzer.Initialize();
            dynamicAnalyzer.Solve();
        }
    }
}
