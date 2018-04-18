using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Problems.Structural.Providers
{
    public class ElementStructuralStiffnessProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(Element element)
        {
            ////return element.K;
            //foreach (var node in element.EmbeddedNodes)
            //{
            //    foreach (var embeddedElement in node.ElementsDictionary.Values)
            //    {
            //        var hostDOFs = element.ElementType.GetDOFTypes(element).SelectMany(d => d).ToArray();
            //        int totalHostDOFs = hostDOFs.Length;
            //        int totalEmbeddedDOFs = embeddedElement.ElementType.GetDOFTypes(element).SelectMany(d => d).Count();
            //        double[] hostShapeFunctions = ((IEmbeddedHostElement)element.ElementType).GetShapeFunctionsForNode(embeddedElement, node);
                    
            //        Matrix2D<double> transformation = new Matrix2D<double>(totalHostDOFs, totalEmbeddedDOFs);
            //        for (int i = 0; i < totalHostDOFs; i++)
            //        {
            //        }
            //    }
            //}

            return element.ElementType.StiffnessMatrix(element);

        }

        #endregion
    }
}
