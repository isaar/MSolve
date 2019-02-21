using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IStaticProvider_v2 : IAnalyzerProvider_v2
    {
        /// <summary>
        /// Returns the submatrix that corresponds to the free freedom degrees of a subdomain. If A is the matrix corresponding 
        /// to all dofs, f denotes free dofs and c denotes constrained dofs then A = [ Aff Acf^T; Acf Acc ]. This method
        /// builds and returns only Aff.
        /// </summary>
        /// <param name="subdomain">The subdomain whose matrix we are interested in.</param>
        IMatrixView CalculateMatrix(ISubdomain_v2 subdomain);

        /// <summary>
        /// Returns the submatrices that corresponds to the free and constrained freedom degrees of a subdomain. If A is the 
        /// matrix corresponding to all dofs, f denotes free dofs and c denotes constrained dofs then A = [ Aff Acf^T; Acf Acc ]. 
        /// This method builds and returns (Aff, Acf, Acc).
        /// </summary>
        /// <param name="subdomain">The subdomain whose matrices we are interested in.</param>
        (IMatrixView matrixFreeFree, IMatrixView matrixConstrFree, IMatrixView matrixConstrConstr) CalculateSubMatrices(
            ISubdomain_v2 subdomain);
    }
}
