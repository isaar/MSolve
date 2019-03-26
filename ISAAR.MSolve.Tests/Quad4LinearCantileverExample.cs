using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Iterative;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public static class Quad4LinearCantileverExample
    {
        [Fact]
        private static void TestQuad4LinearCantileverExample()
        {
            Numerical.LinearAlgebra.VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 3.76;
            double poissonRatio = 0.3779;
            double thickness = 1.0;
            double nodalLoad = 500.0;

            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 10.0, Y = 0.0, Z = 0.0 };
            Node node3 = new Node { ID = 3, X = 10.0, Y = 10.0, Z = 0.0 };
            Node node4 = new Node { ID = 4, X = 0.0, Y = 10.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);
            nodes.Add(node4);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[4].Constraints.Add(new Constraint { DOF = DOFType.Y });

            // Create Quad4 element
            var quad = new Quad4(material) { Thickness = thickness };
            var element = new Element()
            {
                ID = 1,
                ElementType = quad,
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);
            element.AddNode(model.NodesDictionary[3]);
            element.AddNode(model.NodesDictionary[4]);

            // Element Stiffness Matrix
            var a = quad.StiffnessMatrix(element);

            // Add quad element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

            // define loads
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.X });
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[3], DOF = DOFType.X });
            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            // Choose child analyzer -> Child: Linear or NewtonRaphsonNonLinearAnalyzer
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            // Choose parent analyzer -> Parent: Static or Dynamic
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            Assert.Equal(253.132375961535, linearSystems[1].Solution[0], 8);
        }

        [Fact]
        private static void TestQuad4LinearCantileverExample_v2()
        {
            // Model & subdomains
            var model = new Model_v2();
            int subdomainID = 0;
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Materials
            double youngModulus = 3.76;
            double poissonRatio = 0.3779;
            double thickness = 1.0;
            double nodalLoad = 500.0;
            var material = new ElasticMaterial2D_v2(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            // Nodes
            var nodes = new Node_v2[]
            {
                new Node_v2 { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 },
                new Node_v2 { ID = 2, X = 10.0, Y = 0.0, Z = 0.0 },
                new Node_v2 { ID = 3, X = 10.0, Y = 10.0, Z = 0.0 },
                new Node_v2 { ID = 4, X = 0.0, Y = 10.0, Z = 0.0 }
            };
            for (int i = 0; i < nodes.Length; ++i) model.NodesDictionary.Add(i, nodes[i]);


            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);

            var elementWrapper = new Element_v2()
            {
                ID = 0,
                ElementType = factory.CreateElement(CellType.Quad4, nodes)
            };
            elementWrapper.AddNodes(nodes);
            model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
            model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

            //var a = quad.StiffnessMatrix(element);

            // Prescribed displacements
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });
            model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
            model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });

            // Nodal loads
            model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[1], DOF = DOFType.X });
            model.Loads.Add(new Load_v2() { Amount = nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.X });

            // Solver
            var pcgBuilder = new PcgAlgorithm.Builder();
            pcgBuilder.ResidualTolerance = 1E-6;
            pcgBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(0.5);
            var solverBuilder = new PcgSolver.Builder();
            solverBuilder.PcgAlgorithm = pcgBuilder.Build();
            PcgSolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural_v2(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            // Request output
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory_v2(new int[] { 0 });

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog_v2 log = (DOFSLog_v2)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Assert.Equal(253.132375961535, log.DOFValues[0], 8);
        }
    }
}
