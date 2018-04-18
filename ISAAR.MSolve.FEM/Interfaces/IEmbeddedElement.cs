using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Embedding;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedElement
    {
        IList<EmbeddedNode> EmbeddedNodes { get; }
        Dictionary<DOFType, int> GetInternalNodalDOFs(Element element, Node node);
        double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues);
    }
}
