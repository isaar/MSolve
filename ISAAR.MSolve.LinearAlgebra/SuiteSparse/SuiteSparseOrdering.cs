using System;
using System.Collections.Generic;
using System.Text;

//TODO: make sure these match the pinvoke method arguments
namespace ISAAR.MSolve.LinearAlgebra.SuiteSparse
{
    /// <summary>
    /// These will be used during factorization.
    /// </summary>
    public enum SuiteSparseOrdering
    {
        /// <summary>
        /// No reordering
        /// </summary>
        Natural = 0,

        /// <summary>
        /// Let SuiteSparse try all its default algorithms, until it finds a good ordering.
        /// </summary>
        TryDefaults = 1,

        /// <summary>
        /// Use Approximate Minimal Degree algorithm.
        /// </summary>
        AMD = 2
    }

}
