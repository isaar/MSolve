using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalBenchmarkProvatidis_11_2
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            Model Model = CreateModel();
            IVectorView solution = SolveModel(Model);
            Assert.True(CompareResults(solution));
        }

        private static bool CompareResults(IVectorView solution)
        {
            var comparer = new ValueComparer(1E-8);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
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
            model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;

            // Nodes
            int numNodes = 9;
            var nodes = new Node[numNodes];
            nodes[0] = new Node { ID = 0, X = 0.0, Y = 0.0 };
            nodes[1] = new Node { ID = 1, X = 1.0, Y = 0.0 };
            nodes[2] = new Node { ID = 2, X = 2.0, Y = 0.0 };
            nodes[3] = new Node { ID = 3, X = 0.0, Y = 1.0 };
            nodes[4] = new Node { ID = 4, X = 1.0, Y = 1.0 };
            nodes[5] = new Node { ID = 5, X = 2.0, Y = 1.0 };
            nodes[6] = new Node { ID = 6, X = 0.0, Y = 2.0 };
            nodes[7] = new Node { ID = 7, X = 1.0, Y = 2.0 };
            nodes[8] = new Node { ID = 8, X = 2.0, Y = 2.0 };

            for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

            // Elements
            int numElements = 4;
            var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, k));
            var elements = new ThermalElement2D[4];
            elements[0] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[0], nodes[1], nodes[4], nodes[3] });
            elements[1] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[1], nodes[2], nodes[5], nodes[4] });
            elements[2] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[3], nodes[4], nodes[7], nodes[6] });
            elements[3] = elementFactory.CreateElement(CellType.Quad4, new Node[] { nodes[4], nodes[5], nodes[8], nodes[7] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[i] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
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

            return model;
        }

        private static IVectorView SolveModel(Model model)
        {
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}
