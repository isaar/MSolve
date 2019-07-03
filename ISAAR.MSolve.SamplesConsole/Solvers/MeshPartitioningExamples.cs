using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.GMSH;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.SamplesConsole.Solvers
{
    // Fillet and holes meshes.
    public static class MeshPartitioningExamples
    {
        private const string directory = @"C:\Users\Serafeim\Desktop\MeshPartition\";

        public static void PartitionMeshes()
        {
            //PartitionFilletMesh(directory + @"Fillet\fillet_1272.msh", directory + @"Fillet\fillet_subdomains_1272.vtk");
            //PartitionFilletMesh(directory + @"Fillet\fillet_21420.msh", directory + @"Fillet\fillet_subdomains_21420.vtk");
            //PartitionHolesMesh(directory + @"Holes\holes_4442.msh", directory + @"Holes\holes_subdomains_4442.vtk");
            PartitionHolesMesh(directory + @"Holes\holes_13738.msh", directory + @"Holes\holes_subdomains_13738.vtk");
        }

        private static void PartitionFilletMesh(string meshPath, string plotPath)
        {
            (XModel model, IMesh2D<XNode, XContinuumElement2D> mesh) = CreateModel(meshPath);

            // Define subdomain boundaries
            double tol = 1E-13;
            var regions = new Dictionary<int, RectangularRegion2D>();
            regions[0] = new RectangularRegion2D(  0.0,  0.0, 100.0,  75.0, tol);
            regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
            regions[1] = new RectangularRegion2D(275.0,  0.0, 375.0,  75.0, tol);
            regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[2] = new RectangularRegion2D(100.0,  0.0, 275.0,  75.0, tol);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
            regions[3] = new RectangularRegion2D(100.0, 75.0, 187.5, 150.0, tol);
            regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);
            regions[4] = new RectangularRegion2D(187.5, 75.0, 275.0, 150.0, tol);
            regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);

            // Partition mesh into subdomains
            var regionsGeneral = new Dictionary<int, IRegion2D>();
            foreach (var pair in regions) regionsGeneral[pair.Key] = pair.Value;
            var partitioner = new GuidedPartioner2D<XNode, XContinuumElement2D>(mesh, regionsGeneral);
            Dictionary<int, List<XContinuumElement2D>> subdomains = partitioner.CreateSubdomains();

            // Add subdomains to model and define node / element / subdomain connectivity
            for (int s = 0; s < regions.Count; ++s) model.Subdomains.Add(s, new XSubdomain(s));
            foreach (var idElementsPair in subdomains)
            {
                foreach (XContinuumElement2D element in idElementsPair.Value)
                {
                    model.Subdomains[idElementsPair.Key].Elements.Add(element);
                }
            }
            model.ConnectDataStructures();

            // Plot the resulting subdomains
            var writer = new MeshPartitionWriter();
            var nodesPerSubdomain = new Dictionary<int, IReadOnlyList<XNode>>(); 
            var elementsPerSubdomain = new Dictionary<int, IReadOnlyList<IXFiniteElement>>();
            foreach (int subdomainID in model.Subdomains.Keys)
            {
                nodesPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Nodes;
                elementsPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Elements;
            }
            writer.WriteSubdomainElements(plotPath, nodesPerSubdomain, elementsPerSubdomain);
        }

        public static void PartitionHolesMesh(string meshPath, string plotPath)
        {
            (XModel model, IMesh2D<XNode, XContinuumElement2D> mesh) = CreateModel(meshPath);

            // Geometric data
            double tol = 1E-13;
            double minX = 0.0, maxX = 20.0, minY = 0.0, maxY = 10.0;
            double boundaryY = 5.0;
            double leftBoundary1X = 3.5, leftBoundary2X = 7.0, leftBoundary3X = 10.5;
            double rightBoundary1X = 9.5, rightBoundary2X = 13.0, rightBoundary3X = 16.5;

            // Define subdomain boundaries
            var regions = new Dictionary<int, RectangularRegion2D>();

            regions[0] = new RectangularRegion2D(minX, minY, leftBoundary1X, boundaryY, tol);
            regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
            regions[0].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[1] = new RectangularRegion2D(leftBoundary1X, minY, leftBoundary2X, boundaryY, tol);
            regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
            regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[1].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[2] = new RectangularRegion2D(leftBoundary2X, minY, leftBoundary3X, boundaryY, tol);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[2].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[3] = new RectangularRegion2D(leftBoundary3X, minY, maxX, boundaryY, tol);
            regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Up);
            regions[3].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);

            regions[4] = new RectangularRegion2D(minX, boundaryY, rightBoundary1X, maxY, tol);
            regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[4].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[5] = new RectangularRegion2D(rightBoundary1X, boundaryY, rightBoundary2X, maxY, tol);
            regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[5].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[6] = new RectangularRegion2D(rightBoundary2X, boundaryY, rightBoundary3X, maxY, tol);
            regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);
            regions[6].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Right);

            regions[7] = new RectangularRegion2D(rightBoundary3X, boundaryY, maxX, maxY, tol);
            regions[7].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Down);
            regions[7].AddBoundaryEdge(RectangularRegion2D.RectangleEdge.Left);

            // Partition mesh into subdomains
            var regionsGeneral = new Dictionary<int, IRegion2D>();
            foreach (var pair in regions) regionsGeneral[pair.Key] = pair.Value;
            var partitioner = new GuidedPartioner2D<XNode, XContinuumElement2D>(mesh, regionsGeneral);
            Dictionary<int, List<XContinuumElement2D>> subdomains = partitioner.CreateSubdomains();

            // Add subdomains to model and define node / element / subdomain connectivity
            for (int s = 0; s < regions.Count; ++s) model.Subdomains.Add(s, new XSubdomain(s));
            foreach (var idElementsPair in subdomains)
            {
                foreach (XContinuumElement2D element in idElementsPair.Value)
                {
                    model.Subdomains[idElementsPair.Key].Elements.Add(element);
                }
            }
            model.ConnectDataStructures();

            // Plot the resulting subdomains
            var writer = new MeshPartitionWriter();
            var nodesPerSubdomain = new Dictionary<int, IReadOnlyList<XNode>>();
            var elementsPerSubdomain = new Dictionary<int, IReadOnlyList<IXFiniteElement>>();
            foreach (int subdomainID in model.Subdomains.Keys)
            {
                nodesPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Nodes;
                elementsPerSubdomain[subdomainID] = model.Subdomains[subdomainID].Elements;
            }
            writer.WriteSubdomainElements(plotPath, nodesPerSubdomain, elementsPerSubdomain);
        }

        private static (XModel model, IMesh2D<XNode, XContinuumElement2D> mesh) CreateModel(string meshPath)
        {
            // Mesh generation
            var reader = new GmshReader<XNode>(meshPath);
            (IReadOnlyList<XNode> nodes, IReadOnlyList<CellConnectivity<XNode>> elementConnectivities) = reader.CreateMesh(
                (id, x, y, z) => new XNode(id, x, y, z));

            // Nodes
            var model = new XModel();
            foreach (XNode node in nodes) model.Nodes.Add(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            var jIntegration =
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(4, 4));

            // Elements
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(2.1E7, 0.3);
            var factory = new XContinuumElement2DFactory(integration, jIntegration, material);
            var cells = new XContinuumElement2D[elementConnectivities.Count];
            for (int e = 0; e < cells.Length; ++e)
            {
                XContinuumElement2D element = factory.CreateElement(e, CellType.Quad4, elementConnectivities[e].Vertices);
                cells[e] = element;
                model.Elements.Add(element);
            }

            // Mesh usable for crack-mesh interaction
            //var boundary = new FilletBoundary();
            IDomain2DBoundary boundary = null;
            model.Boundary = boundary;
            var mesh = new BidirectionalMesh2D<XNode, XContinuumElement2D>(model.Nodes, cells, boundary);

            return (model, mesh);
        }

        
    }
}
