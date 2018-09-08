using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(IElement element)
        {
            return element.IElementType.MassMatrix(element);
        }
        #endregion
    }
}
