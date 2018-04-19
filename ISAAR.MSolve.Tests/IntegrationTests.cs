using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class IntegrationTests
    {
        [Fact]
        public void TestSolveHexaCantileverBeam()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            HexaSimpleCantileverBeam.MakeCantileverBeam(model, 0, 0, 0, model.NodesDictionary.Count + 1, model.ElementsDictionary.Count + 1, 1);

            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[16], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[17], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[18], DOF = DOFType.Z });
            model.Loads.Add(new Load() { Amount = -0.25, Node = model.Nodes[19], DOF = DOFType.Z });

            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            //NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            double[] expectedDisplacements = new double[]
            {
                -0.0000025899520106, -0.0000004898560318, -0.0000031099520106, -0.0000025899520106, 0.0000004898560318,
                -0.0000031099520106, 0.0000025899520106, 0.0000004898560318, -0.0000031099520106, 0.0000025899520106,
                -0.0000004898560318, -0.0000031099520106, -0.0000045673419128, -0.0000002423136749, -0.0000107872459340,
                -0.0000045673419128, 0.0000002423136749, -0.0000107872459340, 0.0000045673419128, 0.0000002423136749,
                -0.0000107872459340, 0.0000045673419128, -0.0000002423136749, -0.0000107872459340, -0.0000057299058132,
                -0.0000001253780263, -0.0000216044936601, -0.0000057299058132, 0.0000001253780263, -0.0000216044936601,
                0.0000057299058132, 0.0000001253780263, -0.0000216044936601, 0.0000057299058132, -0.0000001253780263,
                -0.0000216044936601, -0.0000061325564473, -0.0000000425738760, -0.0000339869559207, -0.0000061325564473,
                0.0000000425738760, -0.0000339869559207, 0.0000061325564473, 0.0000000425738760, -0.0000339869559207,
                0.0000061325564473, -0.0000000425738760, -0.0000339869559207
            };

            for (int i = 0; i < expectedDisplacements.Length; i++)
                Assert.Equal(expectedDisplacements[i], linearSystems[1].Solution[i], 10);
        }

        [Fact]
        public void SolveCantileverBeam2D()
        {
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            double nodalLoad = 10.0;

            ElasticMaterial material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 5.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);

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
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });


            // Create a new Beam2D element
            var beam = new EulerBeam2D(youngModulus)
            {
                SectionArea = 1,
                MomentOfInertia = .1
            };

            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            var a = beam.StiffnessMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(-2.08333333333333333e-5, linearSystems[1].Solution[1], 10);
        }


        [Fact]
        private static void SolveRandomVariableBeam2DWithMonteCarlo()
        {
            #region Beam2D Geometry Data
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            double nodalLoad = 10.0;

            var coefficientProvider = new RandomVariableTargetEvaluator(1 / youngModulus, 0.1 / youngModulus, RandomVariableDistributionType.Normal);
            StochasticElasticMaterial material = new StochasticElasticMaterial(coefficientProvider)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
            Node node2 = new Node { ID = 2, X = 5.0, Y = 0.0, Z = 0.0 };
            nodes.Add(node1);
            nodes.Add(node2);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
                model.NodesDictionary.Add(i + 1, nodes[i]);

            // Constrain bottom nodes of the model
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = DOFType.RotZ });


            // Create a new Beam2D element
            var beam = new EulerBeam2D(youngModulus)
            {
                SectionArea = 1,
                MomentOfInertia = .1
            };

            var element = new Element()
            {
                ID = 1,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[1]);
            element.AddNode(model.NodesDictionary[2]);

            var a = beam.StiffnessMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();
            #endregion

            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            Analyzers.LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            MonteCarloAnalyzerWithStochasticMaterial stohasticAnalyzer =
                new MonteCarloAnalyzerWithStochasticMaterial(model, provider, parentAnalyzer, linearSystems,
                    coefficientProvider, 1, 100000);
            stohasticAnalyzer.Initialize();
            stohasticAnalyzer.Solve();

            Assert.Equal(-2.08333333333333333e-5, stohasticAnalyzer.MonteCarloMeanValue, 8);
        }

    }
}
