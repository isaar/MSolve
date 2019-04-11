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
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Iterative;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public static class Quad4LinearCantileverExample
    {
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
