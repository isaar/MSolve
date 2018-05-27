using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Mesh.Providers;
using ISAAR.MSolve.XFEM.Geometry.Shapes;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Output.VTK;

namespace ISAAR.MSolve.XFEM.Tests.Subdomains
{
    class AutomaticDecompositionTest
    {
        private const double L = 12.0;
        private const double h = 4.0;

        public static void Run()
        {
            BiMesh2D mesh = CreateMesh();
            PolygonalRegion[] regions = AssignRegions1();
            var decomposer = new GuideDecomposer(regions, mesh);
            XCluster2D cluster = decomposer.CreateSubdomains();
            WriteDecomposition(mesh, regions, cluster);
        }

        private static PolygonalRegion[] AssignRegions1()
        {
            var regions = new PolygonalRegion[4];

            var vertices1 = new CartesianPoint2D[4];
            vertices1[0] = new CartesianPoint2D(0.0, h);
            vertices1[1] = new CartesianPoint2D(0.0, 0.0);
            vertices1[2] = new CartesianPoint2D(L / 5.0, 0.0);
            vertices1[3] = new CartesianPoint2D(L / 3.0, h);
            var boundaries1 = new HashSet<LineSegment2D>();
            boundaries1.Add(new LineSegment2D(vertices1[2], vertices1[3]));
            regions[0] = new PolygonalRegion(vertices1, boundaries1);

            var vertices2 = new CartesianPoint2D[4];
            vertices2[0] = new CartesianPoint2D(L / 3.0, h);
            vertices2[1] = new CartesianPoint2D(L / 5.0, 0.0);
            vertices2[2] = new CartesianPoint2D(L /2.0, 0.0);
            vertices2[3] = new CartesianPoint2D(L / 2.0, h);
            var boundaries2 = new HashSet<LineSegment2D>();
            boundaries2.Add(new LineSegment2D(vertices2[0], vertices2[1]));
            boundaries2.Add(new LineSegment2D(vertices2[2], vertices2[3]));
            regions[1] = new PolygonalRegion(vertices2, boundaries2);

            var vertices3 = new CartesianPoint2D[4];
            vertices3[0] = new CartesianPoint2D(L / 2.0, h);
            vertices3[1] = new CartesianPoint2D(L / 2.0, 0.0);
            vertices3[2] = new CartesianPoint2D(3.0 * L / 4.0, 0.0);
            vertices3[3] = new CartesianPoint2D(2.0 * L / 3.0, h);
            var boundaries3 = new HashSet<LineSegment2D>();
            boundaries3.Add(new LineSegment2D(vertices3[0], vertices3[1]));
            boundaries3.Add(new LineSegment2D(vertices3[2], vertices3[3]));
            regions[2] = new PolygonalRegion(vertices3, boundaries3);

            var vertices4 = new CartesianPoint2D[4];
            vertices4[0] = new CartesianPoint2D(2.0 * L / 3.0, h);
            vertices4[1] = new CartesianPoint2D(3.0 * L / 4.0, 0.0);
            vertices4[2] = new CartesianPoint2D(L, 0.0);
            vertices4[3] = new CartesianPoint2D(L, h);
            var boundaries4 = new HashSet<LineSegment2D>();
            boundaries4.Add(new LineSegment2D(vertices4[0], vertices4[1]));
            regions[3] = new PolygonalRegion(vertices4, boundaries4);

            return regions;
        }

        private static PolygonalRegion[] AssignRegions2()
        {
            var regions = new PolygonalRegion[4];

            var vertices1 = new CartesianPoint2D[4];
            vertices1[0] = new CartesianPoint2D(0.0, 0.0);
            vertices1[1] = new CartesianPoint2D(L / 2.0, 0.0);
            vertices1[2] = new CartesianPoint2D(L / 2.0, h/2);
            vertices1[3] = new CartesianPoint2D(0.0, h/2);
            var boundaries1 = new HashSet<LineSegment2D>();
            boundaries1.Add(new LineSegment2D(vertices1[1], vertices1[2]));
            boundaries1.Add(new LineSegment2D(vertices1[2], vertices1[3]));
            regions[0] = new PolygonalRegion(vertices1, boundaries1);

            var vertices2 = new CartesianPoint2D[4];
            vertices2[0] = new CartesianPoint2D(L / 2.0, 0);
            vertices2[1] = new CartesianPoint2D(L, 0.0);
            vertices2[2] = new CartesianPoint2D(L, h / 2.0);
            vertices2[3] = new CartesianPoint2D(L / 2.0, h / 2.0);
            var boundaries2 = new HashSet<LineSegment2D>();
            boundaries2.Add(new LineSegment2D(vertices2[2], vertices2[3]));
            boundaries2.Add(new LineSegment2D(vertices2[3], vertices2[0]));
            regions[1] = new PolygonalRegion(vertices2, boundaries2);

            var vertices3 = new CartesianPoint2D[4];
            vertices3[0] = new CartesianPoint2D(L / 2.0, h / 2.0);
            vertices3[1] = new CartesianPoint2D(L, h / 2.0);
            vertices3[2] = new CartesianPoint2D(L, h);
            vertices3[3] = new CartesianPoint2D(L / 2.0, h);
            var boundaries3 = new HashSet<LineSegment2D>();
            boundaries3.Add(new LineSegment2D(vertices3[0], vertices3[1]));
            boundaries3.Add(new LineSegment2D(vertices3[3], vertices3[0]));
            regions[2] = new PolygonalRegion(vertices3, boundaries3);

            var vertices4 = new CartesianPoint2D[4];
            vertices4[0] = new CartesianPoint2D(0.0, h / 2.0);
            vertices4[1] = new CartesianPoint2D(L / 2.0, h / 2.0);
            vertices4[2] = new CartesianPoint2D(L / 2.0, h);
            vertices4[3] = new CartesianPoint2D(0.0, h);
            var boundaries4 = new HashSet<LineSegment2D>();
            boundaries4.Add(new LineSegment2D(vertices4[0], vertices4[1]));
            boundaries4.Add(new LineSegment2D(vertices4[1], vertices4[2]));
            regions[3] = new PolygonalRegion(vertices4, boundaries4);

            return regions;
        }

        private static BiMesh2D CreateMesh()
        {
            // Mesh generation
            double elementSize = 0.1;
            var meshGen = new UniformRectilinearMeshGenerator(L, h, (int)(L / elementSize) + 1, (int)(h / elementSize) + 1);
            (XNode2D[] nodes, List<XNode2D[]> elementConnectivity) = meshGen.CreateMesh();

            // Nodes
            var model = new Model2D();
            foreach (XNode2D node in nodes) model.AddNode(node);

            // Integration rules
            var integration = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            var jIntegration = new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order4x4);

            // Elements
            foreach (XNode2D[] elementNodes in elementConnectivity)
            {
                var materialField = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(2e7, 0.3);
                model.AddElement(new XContinuumElement2D(IsoparametricElementType2D.Quad4, elementNodes, materialField,
                    integration, jIntegration));
            }

            // Mesh usable for crack-mesh interaction
            var boundary = new RectangularBoundary(0.0, L, 0.0, h);
            return new BiMesh2D(model.Nodes, model.Elements, boundary);
        }

        private static void WriteDecomposition(BiMesh2D mesh, PolygonalRegion[] regions, XCluster2D cluster)
        {
            string directory = @"C:\Users\Serafeim\Desktop\GRACM\Subdomains\";
            var writer = new DomainDecompositionWriter();

            writer.WriteRegions(directory + "regions.vtk", regions);
            writer.WriteSubdomainElements(directory + "subdomains.vtk", cluster.Subdomains);
            writer.WriteBoundaryNodes(directory + "boundaryNodes", cluster.Subdomains);
        }
    }
}
