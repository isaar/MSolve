using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should this be in the Problems project?
namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementStructuralStiffnessProvider_v2 : IElementMatrixProvider_v2
    {
        public IMatrix Matrix(IElement_v2 element) => element.ElementType.StiffnessMatrix(element);
    }
}
