using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.DomainDecomposition.FETI;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition
{
    public static class Feti1Tests
    {
        [Fact]
        internal static void TestFeti1Solver()
        {
            Model_v2 model = CreateModel();
            Vector globalU = SolveModel(model);
        }

        private static Model_v2 CreateModel()
        {
            // 6 ----- 7 ----- 8
            // |  (3)  |  (4)  |
            // |       |       |
            // 3 ----- 4 ----- 5
            // |  (1)  |  (2)  |
            // |       |       |
            // 0 ----- 1 ----- 2

            // Material
            double thickness = 1.0;
            var material = new ElasticMaterial2D_v2(StressState2D.PlaneStress) { YoungModulus = 2.1E7, PoissonRatio = 0.3 };
            var dynamicProperties = new DynamicMaterial(1.0, 0.0, 0.0);

            // Model with 4 subdomain
            var model = new Model_v2();
            model.SubdomainsDictionary.Add(0, new Subdomain_v2(0));
            model.SubdomainsDictionary.Add(1, new Subdomain_v2(1));
            model.SubdomainsDictionary.Add(2, new Subdomain_v2(2));
            model.SubdomainsDictionary.Add(3, new Subdomain_v2(3));

            // Generate mesh
            double minX = 0.0, minY = 0.0, maxX = 2.0, maxY = 2.0;
            int numElementsX = 2, numElementsY = 2;
            var meshGenerator = new UniformMeshGenerator2D_v2(minX, minY, maxX, maxY, numElementsX, numElementsY);
            (IReadOnlyList<Node_v2> vertices, IReadOnlyList<CellConnectivity_v2> cells) = meshGenerator.CreateMesh();

            // Add nodes to the model
            for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

            // Elements
            double horizontalBoundaryY = minY + (maxY - minY) / 2.0;
            double verticalBoundaryX = minX + (maxX - minX) / 2.0;
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
            for (int e = 0; e < cells.Count; ++e)
            {
                // Domain decomposition
                int subdomainID;
                if (cells[e].Vertices.All(node => (node.X <= verticalBoundaryX) && (node.Y <= horizontalBoundaryY)))
                {
                    subdomainID = 0;
                }
                else if (cells[e].Vertices.All(node => (node.X >= verticalBoundaryX) && (node.Y <= horizontalBoundaryY)))
                {
                    subdomainID = 1;
                }
                else if (cells[e].Vertices.All(node => (node.X <= verticalBoundaryX) && (node.Y >= horizontalBoundaryY)))
                {
                    subdomainID = 2;
                }
                else if (cells[e].Vertices.All(node => (node.X >= verticalBoundaryX) && (node.Y >= horizontalBoundaryY)))
                {
                    subdomainID = 3;
                }
                else throw new Exception("This element does not belong to any subdomain");

                // Create elements
                ContinuumElement2D element = factory.CreateElement(cells[e].CellType, cells[e].Vertices);
                var elementWrapper = new Element_v2() { ID = e, ElementType = element };
                foreach (Node_v2 node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(e, elementWrapper);
                model.SubdomainsDictionary[subdomainID++].Elements.Add(elementWrapper);
            }

            // Boundary conditions
            double tol = 1E-10;
            Node_v2 bottomLeft = model.Nodes.Where(
                node => (Math.Abs(node.X - minX) <= tol) && (Math.Abs(node.Y - minY) <= tol)).First();
            bottomLeft.Constraints.Add(new Constraint() { DOF = DOFType.X, Amount = 0.0 });
            bottomLeft.Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });

            Node_v2 bottomRight = model.Nodes.Where(
                node => (Math.Abs(node.X - maxX) <= tol) && (Math.Abs(node.Y - minY) <= tol)).First();
            bottomRight.Constraints.Add(new Constraint() { DOF = DOFType.Y, Amount = 0.0 });

            return model;
        }

        private static Vector SolveModel(Model_v2 model)
        {
            // Solver
            var solver = new FetiLvl1Solver(model, 1E-7);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.GatherGlobalDisplacements();
        }
    }
}
