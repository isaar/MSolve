using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Providers
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
