using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement
    {
        int ID { get; set; }
        IElementType ElementType { get; }
        IList<INode> Nodes { get; }
        ISubdomain Subdomain { get; }
    }
}
