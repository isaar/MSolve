using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPorousFiniteElement : IFiniteElement
    {
        IMatrix2D PermeabilityMatrix(IElement element);
        IMatrix2D CouplingMatrix(IElement element);
        IMatrix2D SaturationMatrix(IElement element);
    }
}
