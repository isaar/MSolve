using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Problems.Structural.Providers
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
