using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Embedding;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedHostElement_v2
    {
        EmbeddedNode_v2 BuildHostElementEmbeddedNode(Element_v2 element, Node_v2 node,
            IEmbeddedDOFInHostTransformationVector_v2 transformationVector);
        double[] GetShapeFunctionsForNode(Element_v2 element, EmbeddedNode_v2 node); 
    }
}
