using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
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
            model.NodesDictionary[1].Constraints.Add(DOFType.X);
            model.NodesDictionary[1].Constraints.Add(DOFType.Y);
            model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);


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


        //[Fact]
        //private static void SolveRandomVariableBeam2DWithMonteCarlo()
        //{
        //    #region Beam2D Geometry Data
        //    VectorExtensions.AssignTotalAffinityCount();
        //    double youngModulus = 2.0e08;
        //    double poissonRatio = 0.3;
        //    double nodalLoad = 10.0;

        //    var coefficientProvider = new RandomVariableTargetEvaluator(1 / youngModulus, 0.1 / youngModulus, RandomVariableDistributionType.Normal);
        //    StochasticElasticMaterial material = new StochasticElasticMaterial(coefficientProvider)
        //    {
        //        YoungModulus = youngModulus,
        //        PoissonRatio = poissonRatio,
        //    };

        //    // Node creation
        //    IList<Node> nodes = new List<Node>();
        //    Node node1 = new Node { ID = 1, X = 0.0, Y = 0.0, Z = 0.0 };
        //    Node node2 = new Node { ID = 2, X = 5.0, Y = 0.0, Z = 0.0 };
        //    nodes.Add(node1);
        //    nodes.Add(node2);

        //    // Model creation
        //    Model model = new Model();

        //    // Add a single subdomain to the model
        //    model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

        //    // Add nodes to the nodes dictonary of the model
        //    for (int i = 0; i < nodes.Count; ++i)
        //        model.NodesDictionary.Add(i + 1, nodes[i]);

        //    // Constrain bottom nodes of the model
        //    model.NodesDictionary[1].Constraints.Add(DOFType.X);
        //    model.NodesDictionary[1].Constraints.Add(DOFType.Y);
        //    model.NodesDictionary[1].Constraints.Add(DOFType.RotZ);


        //    // Create a new Beam2D element
        //    var beam = new EulerBeam2D(youngModulus)
        //    {
        //        SectionArea = 1,
        //        MomentOfInertia = .1
        //    };

        //    var element = new Element()
        //    {
        //        ID = 1,
        //        ElementType = beam
        //    };

        //    // Add nodes to the created element
        //    element.AddNode(model.NodesDictionary[1]);
        //    element.AddNode(model.NodesDictionary[2]);

        //    var a = beam.StiffnessMatrix(element);

        //    // Add Hexa element to the element and subdomains dictionary of the model
        //    model.ElementsDictionary.Add(element.ID, element);
        //    model.SubdomainsDictionary[1].ElementsDictionary.Add(element.ID, element);

        //    // Add nodal load values at the top nodes of the model
        //    model.Loads.Add(new Load() { Amount = -nodalLoad, Node = model.NodesDictionary[2], DOF = DOFType.Y });

        //    // Needed in order to make all the required data structures
        //    model.ConnectDataStructures();
        //    #endregion

        //    var linearSystems = new Dictionary<int, ILinearSystem>();
        //    linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
        //    SolverSkyline solver = new SolverSkyline(linearSystems[1]);
        //    ProblemStructural provider = new ProblemStructural(model, linearSystems);
        //    Analyzers.LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
        //    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
        //    MonteCarloAnalyzerWithStochasticMaterial stohasticAnalyzer =
        //        new MonteCarloAnalyzerWithStochasticMaterial(model, provider, parentAnalyzer, linearSystems,
        //            coefficientProvider, 1, 100000);
        //    stohasticAnalyzer.Initialize();
        //    stohasticAnalyzer.Solve();

        [Fact]
        private static void SolveQuadCantileverDecompositionTest1()
        {
            #region Quad Cantilever Model
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 3.0e07;
            double poissonRatio = 0.3;
            double nodalLoad = 1000;

            // Create a new elastic 2D material
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            int indexNode = 0;
            for (int i = 0; i < 25; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    model.NodesDictionary.Add(indexNode, new Node()
                    {
                        ID = indexNode++,
                        X = i,
                        Y = j,
                        Z = 0.0
                    });
                }
            }

            // Constrain left nodes of the model
            for (int i = 0; i < 5; i++)
            {
                model.NodesDictionary[i].Constraints.Add(DOFType.X);
                model.NodesDictionary[i].Constraints.Add(DOFType.Y);
            }

            int indexElement = 0;
            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    var element = new Element()
                    {
                        ID = indexElement,
                        ElementType = new Quad4(material)
                    };
                    element.AddNode(model.NodesDictionary[i * 5 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 5 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 5 + j + 1]);
                    element.AddNode(model.NodesDictionary[i * 5 + j + 1]);
                    model.ElementsDictionary.Add(indexElement, element);
                    model.SubdomainsDictionary[0].ElementsDictionary.Add(indexElement++, element);

                }
            }
            // Add nodal load values to node 3
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[124], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            #endregion
            
            AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 3);
            domainDecomposer.UpdateModel();


            Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
            {
                { 0, new int[] {0,4,1,5,8,9,2,6,10,12,13,14,3,7,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}},
                { 1, new int[] {32,36,33,37,40,41,34,38,42,44,45,46,35,39,43,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63}},
                { 2, new int[] {64,68,65,69,72,73,66,70,74,76,77,78,67,71,75,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95}}
            };

            for (int i = 0; i < expectedSubdomains.Count; i++)
            {
                var subdomainElements = model.SubdomainsDictionary[i].ElementsDictionary.Values.ToList();
                Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].ElementsDictionary.Count);
                for (int j = 0; j < expectedSubdomains[i].Length; j++)
                {
                    Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
                }
            }
        }

        [Fact]
        private static void SolveQuadCantileverDecompositionTest2()
        {
            #region Quad Cantilever Model
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 3.0e07;
            double poissonRatio = 0.3;
            double nodalLoad = 1000;

            // Create a new elastic 2D material
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            int indexNode = 0;
            for (int i = 0; i < 25; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    model.NodesDictionary.Add(indexNode, new Node()
                    {
                        ID = indexNode++,
                        X = i,
                        Y = j,
                        Z = 0.0
                    });
                }
            }

            // Constrain left nodes of the model
            for (int i = 0; i < 5; i++)
            {
                model.NodesDictionary[i].Constraints.Add(DOFType.X);
                model.NodesDictionary[i].Constraints.Add(DOFType.Y);
            }

            int indexElement = 0;
            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    var element = new Element()
                    {
                        ID = indexElement,
                        ElementType = new Quad4(material)
                    };
                    element.AddNode(model.NodesDictionary[i * 5 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 5 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 5 + j + 1]);
                    element.AddNode(model.NodesDictionary[i * 5 + j + 1]);
                    model.ElementsDictionary.Add(indexElement, element);
                    model.SubdomainsDictionary[0].ElementsDictionary.Add(indexElement++, element);

                }
            }
            // Add nodal load values to node 3
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[124], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            #endregion

            AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 8);
            domainDecomposer.UpdateModel();


            Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
            {
                { 0, new int[] {0,4,1,5,8,9,2,6,10,12,13,14}},
                { 1, new int[] {3,7,11,15,18,19,17,21,22,23,16,20}},
                { 2, new int[] {24,28,25,29,32,33,26,30,34,36,37,38}},
                { 3, new int[] {27,31,35,39,42,43,41,45,46,47,40,44}},
                { 4, new int[] {48,52,49,53,56,57,50,54,58,60,61,62}},
                { 5, new int[] {51,55,59,63,66,67,65,69,70,71,64,68}},
                { 6, new int[] {72,76,73,77,80,81,74,78,82,84,85,86}},
                { 7, new int[] {75,79,83,87,90,91,89,93,94,95,88,92}}
            };

            for (int i = 0; i < expectedSubdomains.Count; i++)
            {
                var subdomainElements = model.SubdomainsDictionary[i].ElementsDictionary.Values.ToList();
                Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].ElementsDictionary.Count);
                for (int j = 0; j < expectedSubdomains[i].Length; j++)
                {
                    Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
                }
            }
        }

        [Fact]
        private static void SolveQuadCantileverDecompositionTest3()
        {
            #region Quad Cantilever Model
            VectorExtensions.AssignTotalAffinityCount();
            double youngModulus = 3.0e07;
            double poissonRatio = 0.3;
            double nodalLoad = 1000;

            // Create a new elastic 2D material
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Add nodes to the nodes dictonary of the model
            int indexNode = 0;
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    model.NodesDictionary.Add(indexNode, new Node()
                    {
                        ID = indexNode++,
                        X = i,
                        Y = j,
                        Z = 0.0
                    });
                }
            }

            // Constrain left nodes of the model
            for (int i = 0; i < 3; i++)
            {
                model.NodesDictionary[i].Constraints.Add(DOFType.X);
                model.NodesDictionary[i].Constraints.Add(DOFType.Y);
            }

            int indexElement = 0;
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    var element = new Element()
                    {
                        ID = indexElement,
                        ElementType = new Quad4(material)
                    };
                    element.AddNode(model.NodesDictionary[i * 3 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 3 + j]);
                    element.AddNode(model.NodesDictionary[(i + 1) * 3 + j + 1]);
                    element.AddNode(model.NodesDictionary[i * 3 + j + 1]);
                    model.ElementsDictionary.Add(indexElement, element);
                    model.SubdomainsDictionary[0].ElementsDictionary.Add(indexElement++, element);

                }
            }
            // Add nodal load values to node 3
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[17], DOF = DOFType.Y });

            // Needed in order to make all the required data structures
            model.ConnectDataStructures();

            #endregion

            AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 3);
            domainDecomposer.UpdateModel();


            Dictionary<int, int[]> expectedSubdomains = new Dictionary<int, int[]>()
            {
                { 0, new int[] {0,2,1,3}},
                { 1, new int[] {4,6,5,7}},
                { 2, new int[] {8,9}}
            };

            for (int i = 0; i < expectedSubdomains.Count; i++)
            {
                var subdomainElements = model.SubdomainsDictionary[i].ElementsDictionary.Values.ToList();
                Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].ElementsDictionary.Count);
                for (int j = 0; j < expectedSubdomains[i].Length; j++)
                {
                    Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
                }
            }
        }


        [Fact]
        private static void SolveLinearTrussExample()
        {
            VectorExtensions.AssignTotalAffinityCount();

            #region CreateGeometry

            IList<Node> nodes = new List<Node>();
            Node node1 = new Node { ID = 1, X = 0, Y = 0 };
            Node node2 = new Node { ID = 2, X = 0, Y = 40 };
            Node node3 = new Node { ID = 3, X = 40, Y = 40 };

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            double youngMod = 10e6;
            double poisson = 0.3;
            double loadX = 500;
            double loadY = 300;
            double sectionArea = 1.5;
            
            Model trussModel = new Model();

            trussModel.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            for (int i = 0; i < nodes.Count; i++)
            {
                trussModel.NodesDictionary.Add(i + 1, nodes[i]);
            }

            trussModel.NodesDictionary[1].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[1].Constraints.Add(DOFType.Y);
            trussModel.NodesDictionary[2].Constraints.Add(DOFType.X);
            trussModel.NodesDictionary[2].Constraints.Add(DOFType.Y);


            var element1 = new Element() { ID = 1, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };
            var element2 = new Element() { ID = 2, ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea } };

            element1.AddNode(trussModel.NodesDictionary[1]);
            element1.AddNode(trussModel.NodesDictionary[3]);

            element2.AddNode(trussModel.NodesDictionary[2]);
            element2.AddNode(trussModel.NodesDictionary[3]);

            trussModel.ElementsDictionary.Add(element1.ID, element1);
            trussModel.ElementsDictionary.Add(element2.ID, element2);

            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element1.ID, element1);
            trussModel.SubdomainsDictionary[0].ElementsDictionary.Add(element2.ID, element2);

            trussModel.Loads.Add(new Load() { Amount = loadX, Node = trussModel.NodesDictionary[3], DOF = DOFType.X });
            trussModel.Loads.Add(new Load() { Amount = loadY, Node = trussModel.NodesDictionary[3], DOF = DOFType.Y });

            trussModel.ConnectDataStructures();
            #endregion

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically
            linearSystems[0] = new SkylineLinearSystem(0, trussModel.SubdomainsDictionary[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            ProblemStructural provider = new ProblemStructural(trussModel, linearSystems);

            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            Assert.Equal(0.00053333333333333336, linearSystems[0].Solution[0], 10);
            Assert.Equal(0.0017294083664636196, linearSystems[0].Solution[1], 10);
        }
    }
}
