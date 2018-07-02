using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IEmbeddedElement
    {
        //IList<MSolve.FEM.Embedding.EmbeddedNode> EmbeddedNodes { get; }
        Dictionary<DOFType, int> GetInternalNodalDOFs(IElement element, INode node);
        double[] GetLocalDOFValues(IElement hostElement, double[] hostDOFValues);
    }
}
