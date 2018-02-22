using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Embedding;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedHostElement
    {
        EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector);
        double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node); 
    }
}
