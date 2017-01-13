using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedHostElement
    {
        EmbeddedNode BuildHostElementEmbeddedNode(Element element, Node node, IEmbeddedDOFInHostTransformationVector transformationVector);
        double[] GetShapeFunctionsForNode(Element element, EmbeddedNode node); 
    }
}
