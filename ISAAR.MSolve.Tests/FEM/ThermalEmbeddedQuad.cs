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
using ISAAR.MSolve.FEM.Embedding;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalEmbeddedQuad
    {
        private const int subdomainID = 0;
        [Fact]
        public static void ThermalEmbeddedElementExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // Choose model
            ThermalEmbeddedExamplesBuilder.ThermalExampleWithEmbedded(model);
            model.ConnectDataStructures();
            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();

            // Skyline Solver
            var solver = new SolverSkyline(linearSystems[subdomainID]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.SubdomainsDictionary[subdomainID].Forces);
            
            var childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            childAnalyzer.EquivalentLoadsAssemblers = new Dictionary<int, IEquivalentLoadsAssembler>()
            {
                { subdomainID, new EquivalentLoadsAssembler(model.Subdomains[subdomainID], new ElementStructuralStiffnessProvider()) }
            };

            var parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { });

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        public static class ThermalEmbeddedExamplesBuilder
        {
            public static void HostElementsBuilder(Model model)
            {
                int numNodes = 4; 
                var nodes = new Node2D[numNodes];
                nodes[0] = new Node2D(0, 0.0, 0.0);
                nodes[1] = new Node2D(1, 1.0, 0.0);
                nodes[2] = new Node2D(2, 1.0, 1.0);
                nodes[3] = new Node2D(3, 0.0, 1.0);
                for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

                // Dirichlet BC
                model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
                model.NodesDictionary[1].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
                model.NodesDictionary[2].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
                model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
                
                //// Neumann BC
                //double q = 0.0;
                //model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[2], DOF = DOFType.Temperature });
                //model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[5], DOF = DOFType.Temperature });
                //model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = DOFType.Temperature });
               
                // Material
                double density = 1.0;
                double k = 1.0;
                double c = 1.0;
                var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, k));
                var elements = new ThermalElement2D[4];

                // Elements
                int numElements = 1;

                elements[1] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[0], nodes[1], nodes[2], nodes[3] });

                for (int i = 0; i < numElements; ++i)
                {
                    var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                    foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                    model.ElementsDictionary[i] = elementWrapper;
                    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(i, elementWrapper);
                }
            }

            public static void EmbeddedElementsBuilder(Model model)
            {
                // Material
                double density = 1.0;
                double k = 1.0;
                double c = 1.0;

                model.NodesDictionary.Add(4, new Node() { ID = 4, X = 0.25, Y = 0.25});
                model.NodesDictionary.Add(5, new Node() { ID = 5, X = 0.50, Y = 0.50});


            }

        }



        public static void ExampleWithEmbedded(Model model)
        {
            HostElementsBuilder(model);
            EmbeddedElementsBuilder(model);
            var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key == 1).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key == 2).Select(kv => kv.Value), true);
        }
    }
}

