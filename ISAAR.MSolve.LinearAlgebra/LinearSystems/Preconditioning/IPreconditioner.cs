using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning
{
    /// <summary>
    /// Represents a matrix M such that inverse(M) is close to inverse(A) and inverse(M) * A has a smaller condition number 
    /// than A, where A is the original matrix of the linear system A*x=b.
    /// </summary>
    public interface IPreconditioner
    {
        /// <summary>
        /// Applies the preconditioner by solving the system: M * v = w during the preconditioning step of an iterative 
        /// algorithm, where M is the preconditioner and the definition of the vectors v, w depends on the iterative algorithm. 
        /// This is equivalent to evaluating: inverse(M) * w.
        /// </summary>
        /// <param name="rhsVector">The right hand side vector of the system M * v = w.</param>
        /// <exception cref="Exceptions.NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="rhsVector"/> is different than the number of rows of this <see cref="IPreconditioner"/>.
        ///     </exception>
        Vector SolveLinearSystem(Vector rhsVector);
    }
}
