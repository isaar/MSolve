using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.SuiteSparse;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Subdomains
{
    class SubdomainTest1
    {
        public const double L = 40.0;

        public static void Run()
        {
            (Model2D model, ISingleCrack crack) = CreateModel();
            XCluster2D cluster = DefinePartition(model).CreateSubdomains();
            cluster.OrderStandardDofs(model);
            var enrichedSubs = cluster.FindEnrichedSubdomains();
            cluster.DofOrderer.OrderSubdomainDofs(enrichedSubs, enrichedSubs, crack);
            cluster.DofOrderer.WriteToConsole();
            PrintElementDofs(model, cluster);
            BuildSignedBooleanMatrices(cluster);
        }

        public static void BuildSignedBooleanMatrices(XCluster2D cluster)
        {
            var assembler = new XClusterMatrixAssembler();

            SignedBooleanMatrix globalB = assembler.BuildGlobalSignedBooleanMatrix(cluster);
            Console.WriteLine("Global signed boolean matrix (all dofs):");
            var writer = new BooleanMatrixWriter(true);
            writer.WriteToConsole(globalB);
            Console.WriteLine();

            Dictionary<XSubdomain2D, SignedBooleanMatrix> subdomainBs = assembler.BuildSubdomainSignedBooleanMatrices(cluster);
            foreach (var subdomainB in subdomainBs)
            {
                Console.WriteLine($"Signed boolean matrix of subdomain {subdomainB.Key.ID}:");
                writer.WriteToConsole(subdomainB.Value);
                Console.WriteLine();
            }
        }

        public static IDomainDecomposer DefinePartition(Model2D model)
        {
            double tol = 1e-6;
            var regions = new RectangularRegion[4];

            regions[0] = new RectangularRegion(0.0, L / 2, 0.0, L / 2, tol);
            regions[0].AddBoundaryEdge(RectangularRegion.RectangleEdge.Right);
            regions[0].AddBoundaryEdge(RectangularRegion.RectangleEdge.Up);

            regions[1] = new RectangularRegion(L / 2, L, 0.0, L / 2, tol);
            regions[1].AddBoundaryEdge(RectangularRegion.RectangleEdge.Left);
            regions[1].AddBoundaryEdge(RectangularRegion.RectangleEdge.Up);

            regions[2] = new RectangularRegion(0.0, L / 2, L / 2, L, tol);
            regions[2].AddBoundaryEdge(RectangularRegion.RectangleEdge.Right);
            regions[2].AddBoundaryEdge(RectangularRegion.RectangleEdge.Down);

            regions[3] = new RectangularRegion(L / 2, L, L / 2, L, tol);
            regions[3].AddBoundaryEdge(RectangularRegion.RectangleEdge.Left);
            regions[3].AddBoundaryEdge(RectangularRegion.RectangleEdge.Down);

            var decomposer = new NonOverlappingRegionDecomposer2D(model.Nodes, model.Elements, regions);
            return decomposer;
        }

        public static (Model2D, ISingleCrack) CreateModel()
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
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, thickness);
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
            var crackVertex0 = new CartesianPoint2D(0.0, 0.625 * L);
            var crackVertex1 = new CartesianPoint2D(0.375 * L, 0.625 * L);
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
            
            return (model, crack);
        }

        public static void PrintElementDofs(Model2D model, XCluster2D cluster)
        {
            Console.WriteLine();
            for (int e = 0; e < model.Elements.Count; ++e)
            {
                Console.WriteLine($"Element {e} dofs:");
                cluster.DofOrderer.MatchElementToGlobalStandardDofsOf(model.Elements[e], 
                    out IReadOnlyDictionary<int, int> standardMap, out IReadOnlyDictionary<int, int> constrainedMap);
                var enrichedMap = cluster.DofOrderer.MatchElementToGlobalEnrichedDofsOf(model.Elements[e]);
                PrintElementDofMap(standardMap, "global standard");
                Console.WriteLine();
                PrintElementDofMap(enrichedMap, "subdomain enriched");
                Console.WriteLine();
                PrintElementDofMap(constrainedMap, "global constrained");
                Console.WriteLine();
            }
        }

        private static void PrintElementDofMap(IReadOnlyDictionary<int, int> elementToLocal, string localName)
        {
            Console.WriteLine($"Element dof -> {localName} dof:");
            foreach (var pair in elementToLocal)
            {
                Console.WriteLine($"       {pair.Key}    ->       {pair.Value} ");
            }
        }
    }
}
