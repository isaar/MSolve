using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;

//TODO: This should probably follow the design of other writer classes.
namespace ISAAR.MSolve.Logging.VTK
{
    public class VtkMeshPartitionWriter
    {
        public const string vtkReaderVersion = "4.1";

        public void WriteRegions(string path, PolygonalRegion2D[] regions)
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

        public void WriteSubdomainElements<TNode, TElement>(string path, Dictionary<int, IReadOnlyList<TNode>> nodesPerSubdomain,
            Dictionary<int, IReadOnlyList<TElement>> elementsPerSubdomain)
            where TNode: INode
            where TElement : ICell<TNode>
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
                int numSubdomains = nodesPerSubdomain.Count;
                var subdomainNodes = new Dictionary<INode, int>[numSubdomains];
                int numNodes = 0;
                foreach (int s in nodesPerSubdomain.Keys) numNodes += nodesPerSubdomain[s].Count;
                writer.WriteLine($"POINTS {numNodes} double");
                int nodeCounter = 0;
                foreach (int s in nodesPerSubdomain.Keys)
                {
                    subdomainNodes[s] = new Dictionary<INode, int>();
                    foreach (INode node in nodesPerSubdomain[s])
                    {
                        writer.WriteLine($"{node.X} {node.Y} 0.0");
                        subdomainNodes[s][node] = nodeCounter++;
                    }
                }
                writer.WriteLine();

                // Element connectivity
                int numElementData = 0;
                int numElements = 0;
                foreach (int s in nodesPerSubdomain.Keys)
                {
                    numElements += elementsPerSubdomain[s].Count;
                    foreach (TElement element in elementsPerSubdomain[s]) numElementData += 1 + element.Nodes.Count;
                }
                writer.WriteLine($"CELLS {numElements} {numElementData}");
                foreach (int s in nodesPerSubdomain.Keys)
                {
                    foreach (TElement element in elementsPerSubdomain[s])
                    {
                        writer.Write(element.Nodes.Count);
                        foreach (INode node in element.Nodes) writer.Write(" " + subdomainNodes[s][node]);
                        writer.WriteLine();
                    }
                }
                writer.WriteLine();

                // Element types
                writer.WriteLine("CELL_TYPES " + numElements);
                foreach (int s in nodesPerSubdomain.Keys)
                {
                    foreach (TElement element in elementsPerSubdomain[s])
                    {
                        writer.WriteLine(VtkCell.CellTypeCodes[element.CellType]);
                    }
                }
                writer.WriteLine();

                // Assign a different scalar to the nodes of each subdomain
                writer.WriteLine("POINT_DATA " + numNodes);
                writer.WriteLine($"SCALARS subdomainID double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                int id = 0;
                foreach (int s in nodesPerSubdomain.Keys)
                {
                    foreach (var node in nodesPerSubdomain[s])
                    {
                        writer.WriteLine(id);
                    }
                    ++id;
                }
                writer.WriteLine();

                //// Flag the nodes as internal or boundary
                //writer.WriteLine($"SCALARS nodeMembership double 1");
                //writer.WriteLine("LOOKUP_TABLE default");
                //foreach (var subdomain in subdomains)
                //{
                //    foreach (var node in subdomain.AllNodes)
                //    {
                //        writer.WriteLine(subdomain.InternalNodes.Contains(node) ? 1.0 : -1.0);
                //    }
                //}
            }
        }
    }
}
