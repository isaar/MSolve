using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Meshes.Custom;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;

namespace ISAAR.MSolve.SamplesConsole.FEM
{
    public class Cantilever2D
    {
        private const double length = 4.0;
        private const double height = 20.0;
        private const double thickness = 0.1;
        private const double youngModulus = 2E6;
        private const double poissonRatio = 0.3;
        private const double maxLoad = 1000.0; // TODO: this should be triangular

        public static void Run()
        {
            //(Node2D[] nodes, Node2D[][] elementConnectivity) = GenerateMeshManually();
            (Node2D[] nodes, Node2D[][] elementConnectivity) = GenerateUniformMesh();
            Model model = CreateModel(nodes, elementConnectivity);
            SolveLinearStatic(model);
        }

        private static Model CreateModel(Node2D[] nodes, Node2D[][] elementConnectivity)
        {
            // Initialize
            int numNodes = nodes.Length;
            int numElements = elementConnectivity.Length;
            VectorExtensions.AssignTotalAffinityCount();

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
                StressState = "plstress"
            };

            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain() { ID = 0 });

            // Nodes
            for (int i = 0; i < numNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < numElements; ++i)
            {
                ContinuumElement2D element = factory.CreateQuad4(elementConnectivity[i]);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(i, elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node2D[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(DOFType.X);
                constrainedNodes[i].Constraints.Add(DOFType.Y);
            }

            // Loads
            Node2D[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = DOFType.X });

            // Finalize
            model.ConnectDataStructures();
            return model;
        }

        private static (Node2D[] nodes, Node2D[][] elementConnectivity) GenerateMeshManually()
        {
            Node2D[] nodes =
            {
                new Node2D(0, 0.0, 0.0),
                new Node2D(1, length, 0.0),
                new Node2D(2, 0.0, 0.25 * height),
                new Node2D(3, length, 0.25 * height),
                new Node2D(4, 0.0, 0.50 * height),
                new Node2D(5, length, 0.50 * height),
                new Node2D(6, 0.0, 0.75 * height),
                new Node2D(7, length, 0.75 * height),
                new Node2D(8, 0.0, height),
                new Node2D(9, length, height)
            };

            Node2D[][] elementConnectivity =
            {
                new Node2D[] { nodes[0], nodes[1], nodes[3], nodes[2]},
                new Node2D[] { nodes[2], nodes[3], nodes[5], nodes[4]},
                new Node2D[] { nodes[4], nodes[5], nodes[7], nodes[6]},
                new Node2D[] { nodes[6], nodes[7], nodes[9], nodes[8]}
            };

            return (nodes, elementConnectivity);
        }

        private static Model GenerateMeshFromGmsh()
        {
            throw new NotImplementedException();
        }

        private static (Node2D[] nodes, Node2D[][] elementConnectivity) GenerateUniformMesh()
        {
            var meshGen = new UniformMeshGenerator(0.0, 0.0, length, height, 4, 20);
            (Node2D[] nodes, Node2D[][] cellConnectivity) = meshGen.CreateMesh();
            return (nodes, cellConnectivity);
        }

        private static void SolveLinearStatic(Model model)
        {
            // Choose linear equation system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[0]);

            // Choose the provider of the problem -> here a structural problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            // Choose parent and child analyzers -> Parent: Static, Child: Linear
            LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            // Logging displacement, strain, and stress fields.
            string outputDirectory = @"C:\Users\Serafeim\Desktop\Presentation\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
