using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IPorousFiniteElement : IFiniteElement
    {
        IMatrix2D PermeabilityMatrix(Element element);
        IMatrix2D CouplingMatrix(Element element);
        IMatrix2D SaturationMatrix(Element element);
    }
}
