using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement
    {
        int ID { get; set; }
        IElementType ElementType { get; }
        IReadOnlyList<INode> Nodes { get; }
        ISubdomain Subdomain { get; }
    }
}
