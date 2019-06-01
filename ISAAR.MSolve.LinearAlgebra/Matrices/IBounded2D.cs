using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Specifies the dimensions of a matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IBounded2D
    {
        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        int NumRows { get; }
    }
}
