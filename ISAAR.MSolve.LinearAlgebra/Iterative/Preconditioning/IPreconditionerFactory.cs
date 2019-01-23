using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: working with the matrix itself is not always convenient, especially when abstracted behind IMatrixView. E.g. Jacobi
//      preconditioner needs to access the diagonal which is more efficient with the DOK. IncompleteCholesky may also be more 
//      effient on other matrix storage formats than the CSR that will be used for multiplications.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning
{
    /// <summary>
    /// Classes implementing this interface are responsible for the creation of <see cref="IPreconditioner"/> instances.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IPreconditionerFactory
    {
        /// <summary>
        /// Initializes and returns an <see cref="IPreconditioner"/> for the provided <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The original matrix, whose preconditioner will be created.</param>
        IPreconditioner CreatePreconditionerFor(IMatrixView matrix);
    }
}
