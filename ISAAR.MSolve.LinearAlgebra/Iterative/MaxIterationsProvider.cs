using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: use static factory methods instead of constructors.
namespace ISAAR.MSolve.LinearAlgebra.Iterative
{
    /// <summary>
    /// This class decides the maximum number of iterations that will be executed by an iterative algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MaxIterationsProvider
    {
        private readonly int maxIterations = int.MinValue;
        private readonly double maxIterationsOverMatrixOrder = double.MinValue;
        private readonly bool useMaxIterations = false;
        private readonly bool useMaxIterationsOverMatrixOrder = false;

        /// <summary>
        /// Initializes a new <see cref="MaxIterationsProvider"/> instance that will use a fixed number of iterations, regardless
        /// of the matrix order.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations. Must be &gt; 0.</param>
        public MaxIterationsProvider(int maxIterations)
        {
            if (maxIterations < 1) throw new ArgumentException($"Max iterations must be > 0, but were {maxIterations}");
            this.useMaxIterations = true;
            this.maxIterations = maxIterations;
        }

        /// <summary>
        /// Initializes a new <see cref="MaxIterationsProvider"/> instance that will use a percentage of the matrix order as the
        /// max number of iterations.
        /// </summary>
        /// <param name="maxIterationsOverMatrixOrder">
        /// The percentage of the matrix order that will be set as max iterations. It will be rounded up. Constraints:
        /// <paramref name="maxIterationsOverMatrixOrder"/> &gt; 0.0.
        /// </param>
        public MaxIterationsProvider(double maxIterationsOverMatrixOrder)
        {
            if (maxIterationsOverMatrixOrder <= 0.0) throw new ArgumentException(
                $"The ratio of max iterations / matrix order must be > 0.0, but was {maxIterationsOverMatrixOrder}");
            this.useMaxIterationsOverMatrixOrder = true;
            this.maxIterationsOverMatrixOrder = maxIterationsOverMatrixOrder;
        }

        /// <summary>
        /// Gets the number of max iterations that will be performed by an iterative algorithm when solving a linear system with
        /// the provided <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The matrix of the linear system. It must be square.</param>
        public int GetMaxIterationsForMatrix(IIndexable2D matrix)
        {
            if (useMaxIterations) return maxIterations;
            else // use maxIterationsOverMatrixOrder
            {
                Preconditions.CheckSquare(matrix);
                if (maxIterationsOverMatrixOrder == 1.0) return matrix.NumRows;
                else return (int)Math.Ceiling(maxIterationsOverMatrixOrder * matrix.NumRows);
            }
        }
    }
}
