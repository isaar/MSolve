using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;


namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix2D Matrix(IElement element);
    }
}
