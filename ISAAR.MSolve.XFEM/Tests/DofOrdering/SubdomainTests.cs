using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

// TODO: needs simpler subdomain and cluster construction. 
// TODO: needs automatic decomposers. Nothing fancy though.
namespace ISAAR.MSolve.XFEM.Tests.DofOrdering
{
    class SubdomainTests
    {
        private const double L = 40.0;

        public static void Run()
        {
            Model2D model = CreateModel();
            XCluster2D cluster = CreateSubdomains(model);
            cluster.OrderDofs(model);
            cluster.DofOrderer.WriteToConsole();
        }

        private static Model2D CreateModel()
        {
            var model = new Model2D();

            // Material
            double E = 2e7;
            double v = 0.3;
            double thickness = 1.0;

            // Mesh generator
            var meshGen = new UniformRectilinearMeshGenerator(L, L, 4, 4);
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = meshGen.CreateMesh();

            // Nodes
            foreach (XNode2D node in nodes) model.AddNode(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            // Elements
            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlainStress(E, v, thickness);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Prescribed displacements 
            var finder = new EntityFinder(model);
            XNode2D bottomLeft = finder.FindNodeWith(0.0, 0.0);
            XNode2D bottomRight = finder.FindNodeWith(L, 0.0);
            model.AddConstraint(bottomLeft, DisplacementDof.X, 0.0);
            model.AddConstraint(bottomLeft, DisplacementDof.Y, 0.0);
            model.AddConstraint(bottomRight, DisplacementDof.X, 0.05);
            model.AddConstraint(bottomRight, DisplacementDof.Y, 0.0);

            // Prescribed loads
            double load = 1e4;
            XNode2D topLeft = finder.FindNodeWith(0.0, L);
            XNode2D topRight = finder.FindNodeWith(L, L);
            model.AddNodalLoad(topLeft, DisplacementDof.Y, load);
            model.AddNodalLoad(topRight, DisplacementDof.Y, load);

            // Crack
            var boundary = new RectangularBoundary(0.0, L, 0.0, L);
            var mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(model.Nodes, model.Elements, boundary);
            var crackVertex0 = new CartesianPoint2D(0.0, 0.75 * L);
            var crackVertex1 = new CartesianPoint2D(0.375 * L, 0.75 * L);
            var crack = new BasicExplicitCrack2D();
            crack.Mesh = mesh;
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);
            crack.InitializeGeometry(crackVertex0, crackVertex1);

            // Enrichment nodes manually (without tips for simplicity)        
            int[] heavisideNodes = { 10, 11, 12, 15, 16, 17 };
            foreach (var n in heavisideNodes)
            {
                var node = model.Nodes[n];
                double[] enrichmentVal = crack.CrackBodyEnrichment.EvaluateFunctionsAt(node);
                node.EnrichmentItems.Add(crack.CrackBodyEnrichment, enrichmentVal);
            }
            
            return model;
        }

        private static XCluster2D CreateSubdomains(Model2D model)
        {
            // Split the domain in 4 quadrants
            //var finder = new EntityFinder(model); //LINQ might be better for this
            double tol = 1e-6;

            // All nodes of each subdomain
            var all0 = new HashSet<XNode2D>(model.Nodes.Where(node => (node.X <= L / 2 + tol) && (node.Y <= L / 2 + tol)));
            var all1 = new HashSet<XNode2D>(model.Nodes.Where(node => (node.X >= L / 2 - tol) && (node.Y <= L / 2 + tol)));
            var all2 = new HashSet<XNode2D>(model.Nodes.Where(node => (node.X <= L / 2 + tol) && (node.Y >= L / 2 - tol)));
            var all3 = new HashSet<XNode2D>(model.Nodes.Where(node => (node.X >= L / 2 - tol) && (node.Y >= L / 2 - tol)));

            // Boundary nodes of each subdomain
            var boundary01 = model.Nodes.Where(node => (Math.Abs(node.X - L / 2) <= tol) && (node.Y <= L / 2 + tol));
            var boundary23 = model.Nodes.Where(node => (Math.Abs(node.X - L / 2) <= tol) && (node.Y >= L / 2 - tol));
            var boundary02 = model.Nodes.Where(node => (node.X <= L / 2 + tol) && (Math.Abs(node.Y - L / 2) <= tol));
            var boundary13 = model.Nodes.Where(node => (node.X >= L / 2 - tol) && (Math.Abs(node.Y - L / 2) <= tol));
            var boundary0 = new HashSet<XNode2D>(boundary01.Union(boundary02));
            var boundary1 = new HashSet<XNode2D>(boundary01.Union(boundary13));
            var boundary2 = new HashSet<XNode2D>(boundary02.Union(boundary23));
            var boundary3 = new HashSet<XNode2D>(boundary13.Union(boundary23));

            // Create the entities
            var cluster = new XCluster2D();

            var subdomain0 = new XSubdomain2D();
            foreach (var node in boundary0) subdomain0.AddBoundaryNode(node);
            foreach (var node in all0.Except(boundary0)) subdomain0.AddInternalNode(node);
            cluster.AddSubdomain(subdomain0);

            var subdomain1 = new XSubdomain2D();
            foreach (var node in boundary1) subdomain1.AddBoundaryNode(node);
            foreach (var node in all1.Except(boundary1)) subdomain1.AddInternalNode(node);
            cluster.AddSubdomain(subdomain1);

            var subdomain2 = new XSubdomain2D();
            foreach (var node in boundary2) subdomain2.AddBoundaryNode(node);
            foreach (var node in all2.Except(boundary2)) subdomain2.AddInternalNode(node);
            cluster.AddSubdomain(subdomain2);

            var subdomain3 = new XSubdomain2D();
            foreach (var node in boundary3) subdomain3.AddBoundaryNode(node);
            foreach (var node in all3.Except(boundary3)) subdomain3.AddInternalNode(node);
            cluster.AddSubdomain(subdomain3);

            // TODO: add elements in the subdomains!

            return cluster;
        }
    }
}
