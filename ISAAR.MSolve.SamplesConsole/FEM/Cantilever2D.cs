using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.VonMisesStress;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Preprocessor.Meshes.GMSH;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;

namespace ISAAR.MSolve.SamplesConsole.FEM
{
    /// <summary>
    /// A 2D cantilever beam modeled with continuum finite elements.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Cantilever2D
    {
        private const double length = 4.0;
        private const double height = 20.0;
        private const double thickness = 0.1;
        private const double youngModulus = 2E6;
        private const double poissonRatio = 0.3;
        private const double maxLoad = 1000.0; // TODO: this should be triangular

        private const string workingDirectory = @"C:\Users\Serafeim\Desktop\Presentation";
        private static readonly string projectDirectory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName + @"\Resources\GMSH";

        public static void Run()
        {
            // Choose one of the mesh files bundled with the project
            //string meshPath = projectDirectory + "\\cantilever_quad4.msh";
            //string meshPath = projectDirectory + "\\cantilever_quad8.msh";
            //string meshPath = projectDirectory + "\\cantilever_quad9.msh";
            //string meshPath = projectDirectory + "\\cantilever_tri3.msh";
            string meshPath = projectDirectory + "\\cantilever_tri6.msh";

            // Or set a path on your machine
            //string meshPath = @"C:\Users\Serafeim\Desktop\Presentation\cantilever.msh";


            (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateMeshFromGmsh(meshPath);
            //(IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateUniformMesh();
            //(IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) = GenerateMeshManually();

            Model model = CreateModel(nodes, elements);
            //PrintMeshOnly(model);
            SolveLinearStatic(model);
        }

        private static Model CreateModel(IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements)
        {
            // Initialize
            int numNodes = nodes.Count;
            int numElements = elements.Count;
            VectorExtensions.AssignTotalAffinityCount();

            // Materials
            ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
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
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
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

        private static (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) GenerateMeshManually()
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

            CellType2D[] cellTypes = { CellType2D.Quad4, CellType2D.Quad4, CellType2D.Quad4, CellType2D.Quad4 };

            CellConnectivity2D[] elements =
            {
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[0], nodes[1], nodes[3], nodes[2]}),
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[2], nodes[3], nodes[5], nodes[4]}),
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[4], nodes[5], nodes[7], nodes[6]}),
                new CellConnectivity2D(CellType2D.Quad4, new Node2D[] { nodes[6], nodes[7], nodes[9], nodes[8]})
            };

            return (nodes, elements);
        }

        private static (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) GenerateMeshFromGmsh(string path)
        {
            using (var reader = new GmshReader2D(path))
            {
                return reader.CreateMesh();
            }
        }

        private static (IReadOnlyList<Node2D> nodes, IReadOnlyList<CellConnectivity2D> elements) GenerateUniformMesh()
        {
            var meshGen = new UniformMeshGenerator(0.0, 0.0, length, height, 4, 20);
            return meshGen.CreateMesh();
        }

        private static void PrintMeshOnly(Model model)
        {
            var mesh = new VtkMesh2D(model);
            using (var writer = new VtkFileWriter(workingDirectory + "\\mesh.vtk"))
            {
                writer.WriteMesh(mesh.Points, mesh.Cells);
            }
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
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true,
                VonMisesStressCalculator = new PlaneStressVonMises()
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
