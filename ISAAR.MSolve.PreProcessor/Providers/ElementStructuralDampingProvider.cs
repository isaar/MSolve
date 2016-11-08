using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Providers
{
    public class ElementStructuralDampingProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(Element element)
        {
            //return element.M;
            return element.ElementType.DampingMatrix(element);
        }

        #endregion
    }
}
