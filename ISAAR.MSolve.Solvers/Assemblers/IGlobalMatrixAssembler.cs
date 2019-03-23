using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: not sure this interface is required
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the matrix of the linear system that will be solved.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix">The type of the matrix that will be solved.</typeparam>
    public interface IGlobalMatrixAssembler<TMatrix>
        where TMatrix : IMatrix
    {
        /// <summary>
        /// Builds the matrix of the linear system that corresponds to a certain subdomain.
        /// </summary>
        /// <param name="dofOrdering">The freedom degree ordering of the subdomain.</param>
        /// <param name="elements">The (finite) elements of the subdomain.</param>
        /// <param name="matrixProvider">Determines the matrix calculated for each element (e.g. stiffness, mass, etc.)</param>
        TMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements,
            IElementMatrixProvider_v2 matrixProvider);

        /// <summary>
        /// Update internal state when the freedom degree ordering is changed (e.g. reordering, XFEM, adaptive FEM). It 
        /// will be called before modifying the current freedom degree ordering.
        /// </summary>
        void HandleDofOrderingWillBeModified();
    }
}
