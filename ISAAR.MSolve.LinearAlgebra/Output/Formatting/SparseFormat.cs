using System.Collections.Generic;

// TODO: While this works for matrices that use arrays, it doesn't work for DOKs. Ideally the sparse format object should 
// know how to print itself
// TODO: I should also have a triangular format class for Symmetric, Lower and Upper classes.
namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Data transfer object for the internal indexing and value arrays of a sparse matrix and their titles during output.
    /// These arrays cannot be mutated through a <see cref="SparseFormat"/> object.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SparseFormat
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="SparseFormat"/> class.
        /// </summary>
        public SparseFormat()
        {
            this.RawIndexArrays = new Dictionary<string, IReadOnlyList<int>>();
        }

        /// <summary>
        /// Indexing arrays of the sparse matrix and their titles.
        /// </summary>
        public Dictionary<string, IReadOnlyList<int>> RawIndexArrays { get; }

        /// <summary>
        /// The value array of the sparse matrix.
        /// </summary>
        public IReadOnlyList<double> RawValuesArray { get; set; }

        /// <summary>
        /// The title of the value array of the sparse matrix.
        /// </summary>
        public string RawValuesTitle { get; set; } = "Values";
    }
}
