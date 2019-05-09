using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.Decomposition;

//TODO: Move this to Logging project with the rest of VTK classes.
namespace ISAAR.MSolve.XFEM.Output.VTK
{
    internal class DomainDecompositionWriter
    {
        public static string vtkReaderVersion = "4.1";
        private static readonly string directory =
            Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Resources\\";

        public void WriteBoundaryNodes(string path, IReadOnlyList<XSubdomain2D_old> subdomains)
        {
            var boundaryNodes = new Dictionary<INode, double>();
            foreach (var subdomain in subdomains)
            {
                foreach (var node in subdomain.BoundaryNodes) boundaryNodes[node] = 0.0;
            }
            using (var writer = new VtkPointWriter(path))
            {
                writer.WriteScalarField("boundary", boundaryNodes);
            }
        }

        public void WriteRegions(string path, PolygonalRegion[] regions)
        {
            var pointIds = new Dictionary<CartesianPoint, int>();
            var polygons = new HashSet<List<int>>();
            foreach (var region in regions)
            {
                foreach (var point in region.Outline)
                {
                    if (!pointIds.ContainsKey(point)) pointIds.Add(point, pointIds.Count);
                }

                var polygon = new List<int>();
                foreach (var point in region.Outline)
                {
                    polygon.Add(pointIds[point]);
                }
                polygons.Add(polygon);
            }
            var pointsSorted = pointIds.OrderBy(entry => entry.Value);

            using (var writer = new StreamWriter(path))
            {
                // Header
                writer.Write("# vtk DataFile Version ");
                writer.WriteLine(vtkReaderVersion);
                writer.WriteLine(path);
                writer.Write("ASCII\n\n");

                // Vertices 
                writer.WriteLine("DATASET POLYDATA");
                writer.WriteLine($"POINTS {pointIds.Count} double");
                foreach (var point in pointsSorted)
                {
                    writer.WriteLine($"{point.Key.X} {point.Key.Y} 0.0");
                }
                writer.WriteLine();

                // Polygons
                int numData = 0;
                foreach (var polygon in polygons) numData += 1 + polygon.Count;
                writer.WriteLine($"POLYGONS {polygons.Count} {numData}");
                foreach (var polygon in polygons)
                {
                    writer.Write(polygon.Count);
                    foreach (var pointId in polygon)
                    {
                        writer.Write(" " + pointId);
                    }
                    writer.WriteLine();
                }
            }
        }

        public void WriteSubdomainElements(string path, IReadOnlyList<XSubdomain2D_old> subdomains)
        {
            using (var writer = new StreamWriter(path))
            {
                // Header
                writer.Write("# vtk DataFile Version ");
                writer.WriteLine(vtkReaderVersion);
                writer.WriteLine(path);
                writer.Write("ASCII\n\n");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");

                // Write all nodes by repeating the boundary ones
                var subdomainNodes = new Dictionary<XNode, int>[subdomains.Count];
                int numNodes = 0;
                foreach (var subdomain in subdomains) numNodes += subdomain.AllNodes.Count;
                writer.WriteLine($"POINTS {numNodes} double");
                int nodeCounter = 0;
                for (int s = 0; s < subdomains.Count; ++s)
                {
                    subdomainNodes[s] = new Dictionary<XNode, int>();
                    foreach (var node in subdomains[s].AllNodes)
                    {
                        writer.WriteLine($"{node.X} {node.Y} 0.0");
                        subdomainNodes[s][node] = nodeCounter++;
                    }
                }
                writer.WriteLine();

                // Element connectivity
                int numElementData = 0;
                int numElements = 0;
                foreach (var subdomain in subdomains)
                {
                    numElements += subdomain.Elements.Count;
                    foreach (var element in subdomain.Elements) numElementData += 1 + element.Nodes.Count;
                }
                writer.WriteLine($"CELLS {numElements} {numElementData}");
                for (int s = 0; s < subdomains.Count; ++s)
                {
                    foreach (var element in subdomains[s].Elements)
                    {
                        writer.Write(element.Nodes.Count);
                        foreach (var node in element.Nodes) writer.Write(" " + subdomainNodes[s][node]);
                        writer.WriteLine();
                    }
                }
                writer.WriteLine();

                // Element types
                writer.WriteLine("CELL_TYPES " + numElements);
                foreach (var subdomain in subdomains)
                {
                    foreach (var element in subdomain.Elements) writer.WriteLine(VtkCell.CellTypeCodes[element.CellType]);
                }
                writer.WriteLine();

                // Assign a different scalar to the nodes of each subdomain
                writer.WriteLine("POINT_DATA " + numNodes);
                writer.WriteLine($"SCALARS subdomainID double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                int id = 0;
                foreach (var subdomain in subdomains)
                {
                    foreach (var node in subdomain.AllNodes)
                    {
                        writer.WriteLine(id);
                    }
                    ++id;
                }
                writer.WriteLine();

                // Flag the nodes as internal or boundary
                writer.WriteLine($"SCALARS nodeMembership double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                foreach (var subdomain in subdomains)
                {
                    foreach (var node in subdomain.AllNodes)
                    {
                        writer.WriteLine(subdomain.InternalNodes.Contains(node) ? 1.0 : -1.0);
                    }
                }
            }
        }
    }
}
