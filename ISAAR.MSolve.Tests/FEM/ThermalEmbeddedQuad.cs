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
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using System.Diagnostics;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ThermalEmbeddedQuad
    {
        private const int subdomainID = 0;
        private const int hostElementsIDStart = 0;
        private const int embeddedElementsIDStart = 1;
        private const double minX = -0.5, minY = -0.5, maxX = 0.5, maxY = 0.5;
        private const double thickness = 1.0;
        private const int numElementsX = 10, numElementsY = 10;
        private static readonly Vector2 temperatureGradient = Vector2.Create(100.0, 0);
        private const double conductivityMatrix = 1.0, conductivityFiber = 1000.0;

        [Fact]
        public static void ThermalEmbeddedElementExample()
        {
            Model_v2 model = new Model_v2();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain_v2(subdomainID));

            // Choose model
            ThermalExampleWithEmbedded(model);

            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal_v2(model, solver);
            var rve = new ThermalSquareRve(model, Vector2.Create(minX, minY), Vector2.Create(maxX, maxY), thickness, 
                temperatureGradient);
            var homogenization = new HomogenizationAnalyzer(model, solver, provider, rve);

            homogenization.Initialize();
            homogenization.Solve();

            IMatrix conductivity = homogenization.EffectiveConstitutiveTensors[subdomainID];
            Debug.WriteLine($"C = [ {conductivity[0, 0]} {conductivity[0, 1]}; {conductivity[1, 0]} {conductivity[1, 1]}");
        }

        private static void AddHostElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double c = 1.0;

            // Generate mesh
            var meshGenerator = new UniformMeshGenerator2D_v2(minX, minY, maxX, maxY, numElementsX, numElementsY);
            (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add the elements to the model
            var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, conductivityMatrix));
            for (int e = 0; e < cells.Count; ++e)
            {
                ThermalElement2D element = elementFactory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e + hostElementsIDStart, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }
        }

        private static void AddEmbeddedElements(Model_v2 model)
        {
            // Material
            double density = 1.0;
            double c = 1.0;
            double crossSectionArea = 0.1;
            var embeddedMaterial = new ThermalMaterial(density, c, conductivityFiber);

            // Nodes
            int numNonEmbeddedNodes = model.NodesDictionary.Count;
            int embeddedNode1 = numNonEmbeddedNodes + 1; // We do not know if the non embedded node IDs start from 0 or 1. This way there are no duplicate IDs, but there may be a gap.
            int embeddedNode2 = numNonEmbeddedNodes + 2;
            model.NodesDictionary.Add(embeddedNode1, new Node_v2() { ID = embeddedNode1, X = minX, Y = minY + 0.25 });
            model.NodesDictionary.Add(embeddedNode2, new Node_v2() { ID = embeddedNode2, X = maxX, Y = minY + 0.25 });

            // Elements
            Node_v2[] startEndNodes = { model.NodesDictionary[embeddedNode1], model.NodesDictionary[embeddedNode2] };
            var elementType = new ThermalRod(startEndNodes, crossSectionArea, embeddedMaterial);
            int numNonEmbeddedElements = model.ElementsDictionary.Count();
            int embeddedElementID = hostElementsIDStart + numNonEmbeddedElements;
            var elementWrapper = new Element_v2() { ID = embeddedElementID, ElementType = elementType };
            foreach (var node in startEndNodes) elementWrapper.AddNode(node);
            model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
            model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);

            // Apply embedding
            var embeddedGrouping = new ThermalEmbeddedGrouping(model,
                model.ElementsDictionary.Where(x => x.Key < numNonEmbeddedElements).Select(kv => kv.Value),
                model.ElementsDictionary.Where(x => x.Key >= numNonEmbeddedElements).Select(kv => kv.Value),
                new ThermalElementTransformationVector());
            embeddedGrouping.ApplyEmbedding();
        }

        private static void ThermalExampleWithEmbedded(Model_v2 model)
        {
            AddHostElements(model);
            AddEmbeddedElements(model);
        }
    }
}

