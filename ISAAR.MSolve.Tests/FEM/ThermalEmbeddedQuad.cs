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
using System.Linq;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalEmbeddedQuad
    {
        private const int subdomainID = 0;
        private const int hostElementsIDStart = 0;
        private const int embeddedElementsIDStart = 1;


        [Fact]
        public static void ThermalEmbeddedElementExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });

            // Choose model
            ThermalExampleWithEmbedded(model);
            model.ConnectDataStructures();

            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.SubdomainsDictionary[subdomainID].Forces);

            // Skyline Solver
            var solver = new SolverSkyline(linearSystems[subdomainID]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

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

        public static void AddHostElements(Model model)
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
            //model.NodesDictionary[2].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 100.0 });
            //model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });

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
            var elements = new ThermalElement2D[1];

            // Elements
            int numElements = 1;

            elements[0] = elementFactory.CreateElement(CellType2D.Quad4, new Node2D[] { nodes[0], nodes[1], nodes[2], nodes[3] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element() { ID = hostElementsIDStart + i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
            }
        }

        public static void AddEmbeddedElements(Model model)
        {
            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;
            double crossSectionArea = 0.1;

            model.NodesDictionary.Add(4, new Node() { ID = 4, X = 0.25, Y = 0.25 });
            model.NodesDictionary.Add(5, new Node() { ID = 5, X = 0.50, Y = 0.50 });
            Node[] startEndNodes = { model.NodesDictionary[4], model.NodesDictionary[5] };
            var embeddedMaterial = new ThermalMaterial(density, c, k);

            var elementType = new ThermalRod(startEndNodes, crossSectionArea, embeddedMaterial);

            var elementWrapper = new Element() { ID = embeddedElementsIDStart, ElementType = elementType };
            foreach (var node in startEndNodes) elementWrapper.AddNode(node);
            model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(elementWrapper.ID, elementWrapper);

        }

        public static void ThermalExampleWithEmbedded(Model model)
        {
            AddHostElements(model);
            AddEmbeddedElements(model);
            var embeddedGrouping = new ThermalEmbeddedGrouping(model, 
                model.ElementsDictionary.Where(x => x.Key == hostElementsIDStart).Select(kv => kv.Value), 
                model.ElementsDictionary.Where(x => x.Key == embeddedElementsIDStart).Select(kv => kv.Value), 
                new ThermalElementTransformationVector());
            embeddedGrouping.ApplyEmbedding();
        }
    }
}

