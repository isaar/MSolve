using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should this be in the Problems project?
namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.MassMatrix(element);
    }
}
