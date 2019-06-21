using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning;

//TODO: This should probably follow the design of other writer classes.
//TODO: It should be moved to Solvers.
namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public class MeshPartitionWriter
    {
        public const string vtkReaderVersion = "4.1";
        private const int shuffleSeed = 314159;
        private bool shuffleSubdomainColors;

        public MeshPartitionWriter(bool shuffleSubdomainColors = false)
        {
            this.shuffleSubdomainColors = shuffleSubdomainColors;
        }

        public void WriteBoundaryNodes(string path, IStructuralModel model)
        {
            using (var writer = new StreamWriter(path))
            {
                // Header
                writer.Write("# vtk DataFile Version ");
                writer.WriteLine(vtkReaderVersion);
                writer.WriteLine(path);
                writer.Write("ASCII\n\n");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");

                // Write all nodes without repeating the boundary ones
                var boundaryNodes = new List<INode>();
                foreach (INode node in model.Nodes)
                {
                    if (node.SubdomainsDictionary.Count > 1) boundaryNodes.Add(node);
                }
                writer.WriteLine($"POINTS {boundaryNodes.Count} double");
                foreach (INode node in boundaryNodes)
                {
                    writer.WriteLine($"{node.X} {node.Y} {node.Z}");
                }
                writer.WriteLine();

                // Differentiate between crosspoint nodes (more than 2 subdomains)
                writer.WriteLine("POINT_DATA " + boundaryNodes.Count);
                writer.WriteLine($"SCALARS nodes_boundary_crosspoint double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                foreach (INode node in boundaryNodes)
                {
                    writer.WriteLine((node.SubdomainsDictionary.Count == 2) ? 0.0 : 1.0);
                }
            }
        }

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

        public void WriteSpecialNodes(string path, string nodeCategory, IReadOnlyCollection<INode> specialNodes)
        {
            using (var writer = new StreamWriter(path))
            {
                // Header
                writer.Write("# vtk DataFile Version ");
                writer.WriteLine(vtkReaderVersion);
                writer.WriteLine(path);
                writer.Write("ASCII\n\n");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");

                // Write all nodes without repeating the boundary ones
                writer.WriteLine($"POINTS {specialNodes.Count} double");
                foreach (INode node in specialNodes)
                {
                    writer.WriteLine($"{node.X} {node.Y} {node.Z}");
                }
                writer.WriteLine();

                // Give a value to these nodes so that it can be plotted
                writer.WriteLine("POINT_DATA " + specialNodes.Count);
                writer.WriteLine($"SCALARS {nodeCategory} double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                foreach (INode node in specialNodes)
                {
                    writer.WriteLine(1.0);
                }
            }
        }

        public void WriteSubdomainElements(string path, IStructuralModel model)
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
                IReadOnlyList<ISubdomain> subdomains = model.Subdomains;
                int numSubdomains = subdomains.Count;
                var subdomainNodes = new Dictionary<ISubdomain, Dictionary<INode, int>>();
                int numNodes = 0;
                foreach (ISubdomain subdomain in subdomains) numNodes += subdomain.Nodes.Count;
                writer.WriteLine($"POINTS {numNodes} double");
                int nodeCounter = 0;
                foreach (ISubdomain subdomain in subdomains)
                {
                    subdomainNodes[subdomain] = new Dictionary<INode, int>();
                    foreach (INode node in subdomain.Nodes)
                    {
                        writer.WriteLine($"{node.X} {node.Y} {node.Z}");
                        subdomainNodes[subdomain][node] = nodeCounter++;
                    }
                }
                writer.WriteLine();

                // Element connectivity
                int numElementData = 0;
                int numElements = 0;
                foreach (ISubdomain subdomain in subdomains)
                {
                    numElements += subdomain.Elements.Count;
                    foreach (IElement element in subdomain.Elements) numElementData += 1 + element.Nodes.Count;
                }
                writer.WriteLine($"CELLS {numElements} {numElementData}");
                foreach (ISubdomain subdomain in subdomains)
                {
                    foreach (IElement element in subdomain.Elements)
                    {
                        writer.Write(element.Nodes.Count);
                        foreach (INode node in element.Nodes) writer.Write(" " + subdomainNodes[subdomain][node]);
                        writer.WriteLine();
                    }
                }
                writer.WriteLine();

                // Element types
                writer.WriteLine("CELL_TYPES " + numElements);
                foreach (ISubdomain subdomain in subdomains)
                {
                    foreach (IElement element in subdomain.Elements)
                    {
                        CellType cellType = element.ElementType.CellType;
                        bool canBePlotted = VtkCell.CellTypeCodes.TryGetValue(cellType, out int vtkCellCode);
                        if (!canBePlotted) throw new ArgumentException(
                            "The provided subdomains contain elements that cannot be plotted");
                        writer.WriteLine(vtkCellCode);
                    }
                }
                writer.WriteLine();

                // Assign a different scalar to the nodes of each subdomain
                writer.WriteLine("POINT_DATA " + numNodes);
                writer.WriteLine($"SCALARS subdomainID double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                int[] subdomainColors = Enumerable.Range(0, subdomains.Count).ToArray();
                if (shuffleSubdomainColors) ShuffleFischerYates(subdomainColors);
                for (int s = 0; s < subdomains.Count; ++s)
                {
                    foreach (var node in subdomains[s].Nodes)
                    {
                        writer.WriteLine(subdomainColors[s]);
                    }
                }
                writer.WriteLine();
            }
        }

        //TODO: This does not require each IElementType to implement CellType, but is much harder to use. 
        //TODO: IElement, IElementType and ICell need refactoring.
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
                        writer.WriteLine($"{node.X} {node.Y} {node.Z}");
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

        private static int[] ShuffleFischerYates(int[] array)
        {
            Random rng = new Random(shuffleSeed);
            int n = array.Length;        // The number of items left to shuffle (loop invariant).
            while (n > 1)
            {
                int k = rng.Next(n);      // 0 <= k < n.
                n--;                      // n is now the last pertinent index;
                int temp = array[n];     // swap array[n] with array[k] (does nothing if k == n).
                array[n] = array[k];
                array[k] = temp;
            }
            return array;
        }
    }
}
