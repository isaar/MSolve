using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedElement
    {
        IList<EmbeddedNode> EmbeddedNodes { get; }
        Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node);
        double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues);
    }
}
