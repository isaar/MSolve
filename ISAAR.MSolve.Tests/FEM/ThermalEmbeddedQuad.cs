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
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Solvers.Direct;

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
            Model_v2 model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Choose model
            ThermalExampleWithEmbedded(model);

            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemStructural_v2(model, solver);

            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            //childAnalyzer.EquivalentLoadsAssemblers = new Dictionary<int, IEquivalentLoadsAssembler>()
            //{
            //    { subdomainID, new EquivalentLoadsAssembler(model.Subdomains[0], new ElementStructuralStiffnessProvider()) }
            //};

            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        public static void AddHostElements(Model_v2 model)
        {
            int numNodes = 4;
            var nodes = new Node_v2[numNodes];
            nodes[0] = new Node_v2 { ID = 0, X = 0.0, Y = 0.0 };
            nodes[1] = new Node_v2 { ID = 1, X = 1.0, Y = 0.0 };
            nodes[2] = new Node_v2 { ID = 2, X = 1.0, Y = 1.0 };
            nodes[3] = new Node_v2 { ID = 3, X = 0.0, Y = 1.0 };
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

            elements[0] = elementFactory.CreateElement(CellType.Quad4, new Node_v2[] { nodes[0], nodes[1], nodes[2], nodes[3] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element_v2() { ID = hostElementsIDStart + i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }
        }

        public static void AddEmbeddedElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;
            double crossSectionArea = 0.1;

            model.NodesDictionary.Add(4, new Node_v2() { ID = 4, X = 0.25, Y = 0.25 });
            model.NodesDictionary.Add(5, new Node_v2() { ID = 5, X = 0.50, Y = 0.50 });
            Node_v2[] startEndNodes = { model.NodesDictionary[4], model.NodesDictionary[5] };
            var embeddedMaterial = new ThermalMaterial(density, c, k);

            var elementType = new ThermalRod(startEndNodes, crossSectionArea, embeddedMaterial);

            var elementWrapper = new Element_v2() { ID = embeddedElementsIDStart, ElementType = elementType };
            foreach (var node in startEndNodes) elementWrapper.AddNode(node);
            model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
            model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

        }

        public static void ThermalExampleWithEmbedded(Model_v2 model)
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

