namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    /// <summary>
    /// Data transfer object for measurements taken during the execution of a matrix ordering algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ReorderingStatistics
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="ReorderingStatistics"/> class, with the provided data.
        /// </summary>
        /// <param name="supFactorizedNumNonZeros">An upper bound for the number of non zero entries in a subsequent L*L^T 
        ///     factorization.</param>
        /// <param name="numMovedDenseRows">The number of dense rows that are moved to the end during the ordering 
        ///     algorithm.</param>
        public ReorderingStatistics(int supFactorizedNumNonZeros, int numMovedDenseRows)
        {
            this.SupFactorizedNumNonZeros = supFactorizedNumNonZeros;
            this.NumMovedDenseRows = numMovedDenseRows;
        }

        /// <summary>
        /// An upper bound for the number of non zero entries in a subsequent L*L^T factorization.
        /// </summary>
        public int SupFactorizedNumNonZeros { get; }

        /// <summary>
        /// The number of dense rows that are moved to the end during the ordering algorithm.
        /// </summary>
        public int NumMovedDenseRows { get; }
    }
}
