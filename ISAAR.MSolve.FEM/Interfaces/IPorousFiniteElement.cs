using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPorousFiniteElement: IFiniteElement
    {
        IMatrix PermeabilityMatrix(IElement element);
        IMatrix CouplingMatrix(IElement element);
        IMatrix SaturationMatrix(IElement element);
    }
}
