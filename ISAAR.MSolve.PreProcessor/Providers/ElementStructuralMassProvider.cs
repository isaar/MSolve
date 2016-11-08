using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Providers
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(Element element)
        {
            //return element.M;
            return element.ElementType.MassMatrix(element);
        }

        #endregion
    }
}
