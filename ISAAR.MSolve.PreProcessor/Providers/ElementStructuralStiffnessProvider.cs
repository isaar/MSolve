using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Providers
{
    public class ElementStructuralStiffnessProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D<double> Matrix(Element element)
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
