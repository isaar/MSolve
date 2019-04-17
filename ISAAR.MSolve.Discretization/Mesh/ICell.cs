using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Mesh
{
    public interface ICell<TNode>
    {
        IReadOnlyList<TNode> Nodes { get; }
    }
}
