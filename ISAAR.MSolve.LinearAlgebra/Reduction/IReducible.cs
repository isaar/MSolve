namespace ISAAR.MSolve.LinearAlgebra.Reduction
{
    /// <summary>
    /// Defines the main operation of a reduction over the entries of a matrix, vector or other collection.
    /// </summary>
    /// <param name="entry">Each entry of the matrix, vector or other collection.</param>
    /// <param name="aggregator">The result of the reduction so far.</param>
    /// <returns></returns>
    public delegate double ProcessEntry(double entry, double aggregator);

    /// <summary>
    /// Defines the effect of the reduction on zero entries of the matrix, vector or other collection. Since they are all equal
    /// to 0, they can be processed simultaneously.
    /// </summary>
    /// <param name="numZeros">The number of zeros in the matrix, vector or other collection.</param>
    /// <param name="aggregator">The result of the reduction so far.</param>
    /// <returns></returns>
    public delegate double ProcessZeros(int numZeros, double aggregator);

    /// <summary>
    /// A method that will be applied to the scalar result of the reduction, after all zero and non-zero entries have been 
    /// processed.
    /// </summary>
    /// <param name="aggregator">The result of the reduction after all entries are processed.</param>
    /// <returns></returns>
    public delegate double Finalize(double aggregator);

    /// <summary>
    /// Can be the target of a reduction, namely an operation over all entries that produces a scalar result. For more see 
    /// https://en.wikipedia.org/wiki/Fold_(higher-order_function).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IReducible
    {
        /// <summary>
        /// Applies the reduction defined by <paramref name="identityValue"/>, <paramref name="processEntry"/>, 
        /// <paramref name="processZeros"/> and <paramref name="finalize"/> to this vector, matrix or similar collection.
        /// </summary>
        /// <param name="identityValue">A value that will not affect the result of the reduction.</param>
        /// <param name="processEntry">A method that will be applied to each entry of this object.</param>
        /// <param name="processZeros">A method that will be applied to all zero entries of this object simultaneously.</param>
        /// <param name="finalize">A method that will be applied to the scalar result, after all zero and non-zero entries have
        ///     been processed.</param>
        double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize);
    }
}