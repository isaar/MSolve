using System.Collections.Generic;

//TODO: This could be replaced by IElement, but having an ICell generic on the type of node is necessary to abstract the mesh 
//      classes. I also prefer the nodes propert being IReadOnlyList instead of IList, but IList provides .IndexOf() which is 
//      used in embedding. That can be implemented with an extension method though and the whole embedding design is horrible 
//      anyway.
namespace ISAAR.MSolve.Discretization.Mesh
{
    public interface ICell<TNode>
    {
        CellType CellType { get; }
        IReadOnlyList<TNode> Nodes { get; }
    }
}
