using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;


namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix2D Matrix(Element element);
    }
}
