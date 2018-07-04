using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.Implicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Output.VTK
{
    class IntersectedMeshOutput
    {
        private const int triangleVtkCode = 5;
        private readonly TrackingExteriorCrackLSM crack;
        private readonly Model2D model;
        private readonly string pathNoExtension;
        private readonly CartesianTriangulator triangulator;

        //private IReadOnlyList<VtkPoint2D> points;
        //private IReadOnlyList<VtkCell2D> cells;

        public IntersectedMeshOutput(Model2D model, TrackingExteriorCrackLSM crack, string pathNoExtension)
        {
            this.crack = crack;
            this.model = model;
            this.pathNoExtension = pathNoExtension;
            this.triangulator = new CartesianTriangulator();
        }

        public void WriteMeshData(int step)
        {
            // TODO: guess initial capacities from previous steps or from the model
            var allPoints = new List<VtkPoint2D>();
            var allCells = new List<VtkCell2D>();
            int pointCounter = 0;

            foreach (XContinuumElement2D element in model.Elements)
            {
                bool isEnriched = IsEnriched(element);
                if (!isEnriched)
                {
                    var cellPoints = new VtkPoint2D[element.Nodes.Count];
                    for (int p = 0; p < cellPoints.Length; ++p)
                    {
                        cellPoints[p] = new VtkPoint2D(pointCounter++, element.Nodes[p]);
                        allPoints.Add(cellPoints[p]);
                    }
                    allCells.Add(new VtkCell2D(VtkCell2D.CellTypeCodes[element.ElementType], cellPoints));
                }
                else
                {
                    SortedSet<ICartesianPoint2D> triangleVertices = crack.FindTriangleVertices(element);
                    IReadOnlyList<TriangleCartesian2D> triangles = triangulator.CreateMesh(triangleVertices);
                    foreach (TriangleCartesian2D triangle in triangles)
                    {
                        var cellPoints = new VtkPoint2D[3];

                        cellPoints[0] = new VtkPoint2D(pointCounter++, new CartesianPoint2D(triangle.Vertices[0].Coordinates));
                        allPoints.Add(cellPoints[0]);
                        cellPoints[1] = new VtkPoint2D(pointCounter++, new CartesianPoint2D(triangle.Vertices[1].Coordinates));
                        allPoints.Add(cellPoints[1]);
                        cellPoints[2] = new VtkPoint2D(pointCounter++, new CartesianPoint2D(triangle.Vertices[2].Coordinates));
                        allPoints.Add(cellPoints[2]);

                        allCells.Add(new VtkCell2D(triangleVtkCode, cellPoints));
                    }
                }
            }

            using (var writer = new VtkFileWriter2D($"{pathNoExtension}_{step}.vtk"))
            {
                writer.WriteMesh(allPoints, allCells);
            }
        }

        //TODO: this should be available to all XFEM classes.
        private bool IsEnriched(XContinuumElement2D element)
        {
            foreach (XNode2D node in element.Nodes)
            {
                if (node.EnrichmentItems.Count != 0)
                {
                    return true;
                }
            }
            return false;
        }
    }
}
