using System;
using System.Collections.Generic;
using System.Text;

//TODO: Add time measurements, flop measurements, etc. 
//TODO: Alternatively I could use a logger (pull observer) in the algorithm
//TODO: Needs a better name
//TODO: Each algorithm/author outputs something different. Once enough have been implemented/ported, find an appropriate design
//      to unify them. 
namespace ISAAR.MSolve.LinearAlgebra.Iterative
{
    /// <summary>
    /// Data Transfer Object that collects information about the execution and convergence when solving a linear system with 
    /// an iterative algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class IterativeStatistics
    {
        /// <summary>
        /// The name of the iterative algorithm. E.g CG, PCG, ...
        /// </summary>
        public string AlgorithmName { get; set; }

        /// <summary>
        /// True if the requested tolerance of the residual vector r = b - A * x has been achieved. False if the algorithm has
        /// terminated due to other criteria.
        /// </summary>
        public bool HasConverged { get; set; }

        /// <summary>
        /// The number of iterations that were run before the algorithm terminated. 
        /// </summary>
        public int NumIterationsRequired { get; set; }

        /// <summary>
        /// The value of norm2(b-A*x) / norm2(b-A*x0), where x is the solution vector after the final iteration of the algorithm
        /// and x0 the intial guess for the same solution vector (usually x0 = 0). Depending on the algorithm estimations of this
        /// may be used to check convergence, e.g. PCG uses norm2(r^T*inv(M)*r)/norm(r0^T*inv(M)*r0). Those methods report the
        /// ratio that is checked.
        /// </summary>
        public double ResidualNormRatioEstimation { get; set; }

        /// <summary>
        /// Reports the accumulated data of this <see cref="CGStatistics"/> instance. 
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(AlgorithmName);
            string converged = HasConverged ? " converged successfully." : " failed to converge.";
            sb.Append(converged);
            sb.Append($" A total of {NumIterationsRequired} iterations were run and");
            sb.Append($" norm2(rhs - matrix * xSolution) / norm2(rhs - matrix * xInit) = {ResidualNormRatioEstimation}.");
            return sb.ToString();
        }
    }
}
