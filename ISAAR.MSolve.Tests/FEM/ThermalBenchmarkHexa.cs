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
    public class ThermalBenchmarkHexa
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
            var comparer = new ValueComparer(1E-3);

            //                                         dofs:       4,       5,       6,       7,       8,       9,      13,      14,      15,      16,      17,      18,      22,      23,      24,      25,      26,  27
            var expectedSolution = new Vector(new double[] { 135.054, 158.824, 135.054, 469.004, 147.059, 159.327, 178.178, 147.299, 139.469, 147.059, 191.717, 147.059, 135.054, 158.824, 135.054, 469.004, 147.059, 159.327 });
            int numFreeDofs = 18;
            if (solution.Length != 18) return false;
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
            int numNodes = 27;
            var nodes = new Node3D[numNodes];
            nodes[0] = new Node3D(0, 2.0, 2.0, 2.0);
            nodes[1] = new Node3D(1, 2.0, 1.0, 2.0);
            nodes[2] = new Node3D(2, 2.0, 0.0, 2.0);
            nodes[3] = new Node3D(3, 2.0, 2.0, 1.0);
            nodes[4] = new Node3D(4, 2.0, 1.0, 1.0);
            nodes[5] = new Node3D(5, 2.0, 0.0, 1.0);
            nodes[6] = new Node3D(6, 2.0, 2.0, 0.0);
            nodes[7] = new Node3D(7, 2.0, 1.0, 0.0);
            nodes[8] = new Node3D(8, 2.0, 0.0, 0.0);
            nodes[9] = new Node3D(9, 1.0, 2.0, 2.0);
            nodes[10] = new Node3D(10, 1.0, 1.0, 2.0);
            nodes[11] = new Node3D(11, 1.0, 0.0, 2.0);
            nodes[12] = new Node3D(12, 1.0, 2.0, 1.0);
            nodes[13] = new Node3D(13, 1.0, 1.0, 1.0);
            nodes[14] = new Node3D(14, 1.0, 0.0, 1.0);
            nodes[15] = new Node3D(15, 1.0, 2.0, 0.0);
            nodes[16] = new Node3D(16, 1.0, 1.0, 0.0);
            nodes[17] = new Node3D(17, 1.0, 0.0, 0.0);
            nodes[18] = new Node3D(18, 0.0, 2.0, 2.0);
            nodes[19] = new Node3D(19, 0.0, 1.0, 2.0);
            nodes[20] = new Node3D(20, 0.0, 0.0, 2.0);
            nodes[21] = new Node3D(21, 0.0, 2.0, 1.0);
            nodes[22] = new Node3D(22, 0.0, 1.0, 1.0);
            nodes[23] = new Node3D(23, 0.0, 0.0, 1.0);
            nodes[24] = new Node3D(24, 0.0, 2.0, 0.0);
            nodes[25] = new Node3D(25, 0.0, 1.0, 0.0);
            nodes[26] = new Node3D(26, 0.0, 0.0, 0.0);

            for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

            // Elements
            int numElements = 8;
            var elementFactory = new ThermalElement3DFactory(new ThermalMaterial(density, c, k));
            var elements = new ThermalElement3D[8];
            elements[0] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[13], nodes[4], nodes[3], nodes[12], nodes[10], nodes[1], nodes[0], nodes[9] });
            elements[1] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[14], nodes[5], nodes[4], nodes[13], nodes[11], nodes[2], nodes[1], nodes[10] });
            elements[2] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[16], nodes[7], nodes[6], nodes[15], nodes[13], nodes[4], nodes[3], nodes[12] });
            elements[3] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[17], nodes[8], nodes[7], nodes[16], nodes[14], nodes[5], nodes[4], nodes[13] });
            elements[4] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[22], nodes[13], nodes[12], nodes[21], nodes[19], nodes[10], nodes[9], nodes[18] });
            elements[5] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[23], nodes[14], nodes[13], nodes[22], nodes[20], nodes[11], nodes[10], nodes[19] });
            elements[6] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[25], nodes[16], nodes[15], nodes[24], nodes[22], nodes[13], nodes[12], nodes[21] });
            elements[7] = elementFactory.CreateElement(CellType3D.Hexa8, new Node3D[] { nodes[26], nodes[17], nodes[16], nodes[25], nodes[23], nodes[14], nodes[13], nodes[22] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[i] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(i, elementWrapper);
            }

            // Dirichlet BC
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[1].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[2].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[9].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[10].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[11].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[18].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[19].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            model.NodesDictionary[20].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });

            // Neumann BC
            double q = 100;
            model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[6], DOF = DOFType.Temperature });
            model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[24], DOF = DOFType.Temperature });

            model.ConnectDataStructures();
            return model;
        }

        private static IVector SolveModel(Model model)
        {
            VectorExtensions.AssignTotalAffinityCount();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.SubdomainsDictionary[subdomainID].Forces);
            var solver = new SolverSkyline(linearSystems[subdomainID]);

            var provider = new ProblemStructural(model, linearSystems);

            var childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            childAnalyzer.EquivalentLoadsAssemblers = new Dictionary<int, IEquivalentLoadsAssembler>()
            {
                { subdomainID, new EquivalentLoadsAssembler(model.Subdomains[0], new ElementStructuralStiffnessProvider()) }
            };

            var parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return linearSystems[subdomainID].Solution;
        }
    }
}
