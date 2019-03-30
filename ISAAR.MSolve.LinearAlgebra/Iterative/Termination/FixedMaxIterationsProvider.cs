using System;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.Termination
{
    /// <summary>
    /// A <see cref="IMaxIterationsProvider"/> implementation that will use a fixed number of iterations, regardless of the 
    /// matrix order.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FixedMaxIterationsProvider : IMaxIterationsProvider
    {
        private readonly int maxIterations;

        /// <summary>
        /// Initializes a new instance of <see cref="FixedMaxIterationsProvider"/> with the specified settings.
        /// </summary>
        /// <param name="maxIterations">The fixed number of maximum allowed iterations of the algorithm. Must be &gt; 0.</param>
        public FixedMaxIterationsProvider(int maxIterations)
        {
            if (maxIterations < 1) throw new ArgumentException($"Max iterations must be > 0, but were {maxIterations}");
            this.maxIterations = maxIterations;
        }

        /// <summary>
        /// See <see cref="IMaxIterationsProvider.GetMaxIterations(int)"/>.
        /// </summary>
        public int GetMaxIterations(int matrixOrder) => maxIterations;
    }
}
