using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.Logging.VTK
{
    /// <summary>
    /// VTK unstructured grid. Usually there is an 1-1 correspondence between this and the finite element mesh used for the 
    /// analysis, but sometimes having a different representation is beneficial (e.g. cells that result from intersections of
    /// the finite elements by a crack).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkMesh2D
    {
        private readonly Model model;
        private readonly VtkPoint2D[] points;
        private readonly VtkCell2D[] cells;
        private readonly Dictionary<Node, VtkPoint2D> nodes2Points;
        private readonly Dictionary<ContinuumElement2D, VtkCell2D> elements2Cells;

        public VtkMesh2D(Model model)
        {
            this.model = model;

            IList<Node> nodes = model.Nodes;
            points = new VtkPoint2D[nodes.Count];
            nodes2Points = new Dictionary<Node, VtkPoint2D>();
            for (int i = 0; i < points.Length; ++i)
            {
                Node node = nodes[i];
                points[i] = new VtkPoint2D(i, node.X, node.Y);
                nodes2Points[node] = points[i];
            }

            IList<Element> elements = model.Elements;
            cells = new VtkCell2D[elements.Count];
            elements2Cells = new Dictionary<ContinuumElement2D, VtkCell2D>();
            for (int i = 0; i < cells.Length; ++i)
            {
                ContinuumElement2D element = (ContinuumElement2D)(elements[i].ElementType);
                bool exists = VtkCell2D.cellTypeCodes.TryGetValue(element.Interpolation.CellType, out int code);
                if (!exists) throw new NotImplementedException("Cannot plot elements of type " + element.Interpolation);
                var vertices = new VtkPoint2D[element.Nodes.Count];
                for (int j = 0; j < vertices.Length; ++j) vertices[j] = nodes2Points[element.Nodes[j]];
                cells[i] = new VtkCell2D(code, vertices);
                elements2Cells[element] = cells[i];
            }
        }

        public IReadOnlyList<VtkPoint2D> Points { get { return points; } }
        public IReadOnlyList<VtkCell2D> Cells { get { return cells; } }
        public IReadOnlyDictionary<Node, VtkPoint2D> Nodes2Points { get { return nodes2Points; } }
        public IReadOnlyDictionary<ContinuumElement2D, VtkCell2D> Elements2Celss { get { return elements2Cells; } }
    }
}
