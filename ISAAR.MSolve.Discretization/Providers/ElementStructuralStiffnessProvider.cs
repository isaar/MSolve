using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should this be in the Problems project?
namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementStructuralStiffnessProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.StiffnessMatrix(element);
    }
}
