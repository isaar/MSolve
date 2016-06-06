using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IEmbeddedHostElement
    {
        EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector);
        double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node); 
    }
}
