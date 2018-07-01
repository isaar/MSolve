using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementStructuralStiffnessProvider : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members

        public IMatrix2D Matrix(IElement element)
        {
            return element.IElementType.StiffnessMatrix(element);
        }

        #endregion
    }
}
