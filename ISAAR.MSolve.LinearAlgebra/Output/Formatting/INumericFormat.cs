//TODO: Perhaps I should also specify culture
namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Specifies formatting and alignment for <see cref="double"/> entries of a matrix or vector.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface INumericFormat
    {
        /// <summary>
        /// Returns a tring that will be used as a template for formatting numbers during output.
        /// </summary>
        string GetRealNumberFormat();
    }
}
