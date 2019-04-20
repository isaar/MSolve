using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;

//TODO: The whole conversion INode, IElement -> VtkPoint, VtkCell should not be necessary. It just doubles the memory 
//      requirement for the mesh. Instead INode, IElement should be used. The conversion is necessary for discontinuous meshes
//      though.
namespace ISAAR.MSolve.Logging.VTK
{
    public class VtkMesh<TNode> : IVtkMesh where TNode:INode
    {
        public VtkMesh(IReadOnlyList<TNode> nodes, IReadOnlyList<ICell<TNode>> elements)
        {
            this.OriginalNodes = nodes;
            this.OriginalElements = elements;

            var vtkPoints = new VtkPoint[nodes.Count];
            var nodes2Points = new Dictionary<TNode, VtkPoint>();
            for (int i = 0; i < vtkPoints.Length; ++i)
            {
                TNode node = nodes[i];
                vtkPoints[i] = new VtkPoint(i, node.X, node.Y, node.Z);
                nodes2Points[node] = vtkPoints[i]; //TODO: Even more memory waste.
            }
            this.VtkPoints = vtkPoints;

            var vtkCells = new VtkCell[elements.Count];
            for (int i = 0; i < vtkCells.Length; ++i)
            {
                ICell<TNode> element = elements[i];
                var vertices = new VtkPoint[element.Nodes.Count];
                for (int j = 0; j < vertices.Length; ++j) vertices[j] = nodes2Points[element.Nodes[j]];
                vtkCells[i] = new VtkCell(element.CellType, vertices);
            }
            this.VtkCells = vtkCells;
        }

        public IReadOnlyList<ICell<TNode>> OriginalElements { get; }
        public IReadOnlyList<TNode> OriginalNodes { get; }
        public IReadOnlyList<VtkCell> VtkCells { get; }
        public IReadOnlyList<VtkPoint> VtkPoints { get; }
    }
}
