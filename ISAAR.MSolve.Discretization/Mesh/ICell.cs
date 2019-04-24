using System.Collections.Generic;

//TODO: This could be replaced by IElement, but having an ICell generic on the type of node is necessary to abstract the mesh 
//      classes.
//TODO: If I manage to remove the generics from this one, then the generics can be removed from a lot of other classes. 
namespace ISAAR.MSolve.Discretization.Mesh
{
    public interface ICell<TNode>
    {
        CellType CellType { get; }
        IReadOnlyList<TNode> Nodes { get; }
    }
}
