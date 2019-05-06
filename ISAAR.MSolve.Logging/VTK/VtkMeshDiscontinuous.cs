using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;

namespace ISAAR.MSolve.Logging.VTK
{
    public class VtkMeshDiscontinuous<TNode> : IVtkMesh where TNode:INode
    {
        public VtkMeshDiscontinuous(IReadOnlyList<TNode> nodes, IReadOnlyList<ICell<TNode>> elements)
        {
            this.OriginalNodes = nodes;
            this.OriginalElements = elements;
            var vtkPoints = new List<VtkPoint>();
            var vtkCells = new VtkCell[elements.Count];
            int pointID = 0;

            for (int e = 0; e < elements.Count; ++e)
            {
                ICell<TNode> element = elements[e];
                var cellVertices = new VtkPoint[element.Nodes.Count];
                for (int i = 0; i < element.Nodes.Count; ++i)
                {
                    var point = new VtkPoint(pointID++, element.Nodes[i]);
                    cellVertices[i] = point;
                    vtkPoints.Add(point);
                }
                var cell = new VtkCell(element.CellType, cellVertices);
                vtkCells[e] = cell;
            }
            this.VtkPoints = vtkPoints;
            this.VtkCells = vtkCells;
        }

        public IReadOnlyList<ICell<TNode>> OriginalElements { get; }
        public IReadOnlyList<TNode> OriginalNodes { get; }
        public IReadOnlyList<VtkCell> VtkCells { get; }
        public IReadOnlyList<VtkPoint> VtkPoints { get; }
    }
}
