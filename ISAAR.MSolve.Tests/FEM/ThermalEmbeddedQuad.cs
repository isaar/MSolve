using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

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
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            // Choose model
            ThermalExampleWithEmbedded(model);

            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal(model, solver);
            var rve = new ThermalSquareRve(model, Vector2.Create(minX, minY), Vector2.Create(maxX, maxY), thickness, 
                temperatureGradient);
            var homogenization = new HomogenizationAnalyzer(model, solver, provider, rve);

            homogenization.Initialize();
            homogenization.Solve();

            IMatrix conductivity = homogenization.EffectiveConstitutiveTensors[subdomainID];
            Debug.WriteLine($"C = [ {conductivity[0, 0]} {conductivity[0, 1]}; {conductivity[1, 0]} {conductivity[1, 1]}");
        }

        private static void AddHostElements(Model model)
        {
            // Material
            double density = 1.0;
            double c = 1.0;

            // Generate mesh
            var meshGenerator = new UniformMeshGenerator2D(minX, minY, maxX, maxY, numElementsX, numElementsY);
            (IReadOnlyList<Node> vertices, IReadOnlyList<CellConnectivity> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Add the elements to the model
            var elementFactory = new ThermalElement2DFactory(1.0, new ThermalMaterial(density, c, conductivityMatrix));
            for (int e = 0; e < cells.Count; ++e)
            {
                ThermalElement2D element = elementFactory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element() { ID = e + hostElementsIDStart, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(elementWrapper.ID, elementWrapper);
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }
        }

        private static void AddEmbeddedElements(Model model)
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
            model.NodesDictionary.Add(embeddedNode1, new Node() { ID = embeddedNode1, X = minX, Y = minY + 0.25 });
            model.NodesDictionary.Add(embeddedNode2, new Node() { ID = embeddedNode2, X = maxX, Y = minY + 0.25 });

            // Elements
            Node[] startEndNodes = { model.NodesDictionary[embeddedNode1], model.NodesDictionary[embeddedNode2] };
            var elementType = new ThermalRod(startEndNodes, crossSectionArea, embeddedMaterial);
            int numNonEmbeddedElements = model.ElementsDictionary.Count();
            int embeddedElementID = hostElementsIDStart + numNonEmbeddedElements;
            var elementWrapper = new Element() { ID = embeddedElementID, ElementType = elementType };
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

        private static void ThermalExampleWithEmbedded(Model model)
        {
            AddHostElements(model);
            AddEmbeddedElements(model);
        }
    }
}

