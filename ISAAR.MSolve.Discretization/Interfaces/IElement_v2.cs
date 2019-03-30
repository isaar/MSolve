using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement_v2
    {
        int ID { get; set; }
        IElementType_v2 ElementType { get; }
        IList<INode> Nodes { get; }
        ISubdomain_v2 Subdomain { get; }
    }
}
