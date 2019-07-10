using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Since this interface defines only a single method and the implementations are trivial, I think this interface should be 
//      integrated into a provider interface.
namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElementMatrixProvider
    {
        IMatrix Matrix(IElement element);
    }
}
