using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix2D Matrix(Element element);
    }
}
