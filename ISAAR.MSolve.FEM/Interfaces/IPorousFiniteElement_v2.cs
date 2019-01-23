using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPorousFiniteElement_v2: IFiniteElement_v2
    {
        IMatrix PermeabilityMatrix(IElement_v2 element);
        IMatrix CouplingMatrix(IElement_v2 element);
        IMatrix SaturationMatrix(IElement_v2 element);
    }
}
