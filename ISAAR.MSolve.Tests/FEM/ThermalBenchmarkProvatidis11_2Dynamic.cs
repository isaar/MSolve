using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;
namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalBenchmarkProvatidis11_2Dynamic
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel();
            IVector solution = SolveModel(model);
            Assert.True(CompareResults(solution));
        }

        private static bool CompareResults(IVector solution)
        {
            var comparer = new ValueComparer(1E-5);

            //                                         dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = new Vector(new double[] { 150, 200, 150, 200, 150, 200 });
            int numFreeDofs = 6;
            if (solution.Length != 6) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }

        private static Model CreateModel()
        {
            var model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = subdomainID });

            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;

            // Nodes
            int numNodes = 9;
            var nodes = new Node2D[numNodes];
            nodes[0] = new Node2D(0, 0.0, 0.0);
            nodes[1] = new Node2D(1, 1.0, 0.0);
            nodes[2] = new Node2D(2, 2.0, 0.0);
            nodes[3] = new Node2D(3, 0.0, 1.0);
            nodes[4] = new Node2D(4, 1.0, 1.0);
            nodes[5] = new Node2D(5, 2.0, 1.0);
            nodes[6] = new Node2D(6, 0.0, 2.0);
            nodes[7] = new Node2D(7, 1.0, 2.0);
            nodes[8] = new Node2D(8, 2.0, 2.0);

            for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

            // Elements
            int numElements = 4;
            var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, k));
            var elements = new ThermalElement2D[4];
            elements[0] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[0], nodes[1], nodes[4], nodes[3] });
            elements[1] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[1], nodes[2], nodes[5], nodes[4] });
            elements[2] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[3], nodes[4], nodes[7], nodes[6] });
            elements[3] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[4], nodes[5], nodes[8], nodes[7] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[i] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(i, elementWrapper);
            }

            // Dirichlet BC
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[6].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });

            // Neumann BC
            double q = 50.0;
            model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[2], DOF = DOFType.Temperature });
            model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[5], DOF = DOFType.Temperature });
            model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = DOFType.Temperature });


            ////////model.TimeDependentNodalLoads.Add(new SteadyNodalLoad(q / 2.0) { Node = model.NodesDictionary[2], DOF = DOFType.Temperature });

            model.ConnectDataStructures();
            return model;
        }

        private static IVector SolveModel(Model model)
        {
            VectorExtensions.AssignTotalAffinityCount();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.SubdomainsDictionary[subdomainID].Forces);
            var solver = new SolverSkyline(linearSystems[subdomainID]);

            var provider = new ProblemThermal(model, linearSystems);

            var childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            childAnalyzer.EquivalentLoadsAssemblers = new Dictionary<int, IEquivalentLoadsAssembler>()
            {
                { subdomainID, new EquivalentLoadsAssembler(model.Subdomains[0], new ElementStructuralStiffnessProvider()) }
            };

            var parentAnalyzer = new ThermalDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.5, 0.5, 1000);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return linearSystems[subdomainID].Solution;
        }
    }
}

