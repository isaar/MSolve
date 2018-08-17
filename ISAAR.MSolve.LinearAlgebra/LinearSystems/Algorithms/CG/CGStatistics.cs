using System.Text;

//TODO: Add time measurements, flop measurements, etc. 
//TODO: Alternatively I could use a logger (pull observer) in the algorithm
//TODO: Needs a better name
//TODO: Each algorithm/author outputs something different. Once enough have been implemented/ported, find an appropriate design
//      to unify them. 
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG
{
    /// <summary>
    /// Data Transfer Object that collects information about the execution and convergence when solving a linear system with 
    /// the Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CGStatistics
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
        public int IterationsRequired { get; set; }

        /// <summary>
        /// The value of norm2(b-A*x) / norm2(b-A*x0), where x is the solution vector after the final iteration of the algorithm
        /// and x0 the intial guess for the same solution vector (usually x0 = 0).
        /// </summary>
        public double NormRatio { get; set; }

        /// <summary>
        /// Reports the accumulated data of this <see cref="CGStatistics"/> instance. 
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(AlgorithmName);
            string converged = HasConverged ? " converged successfully." : " failed to converge.";
            sb.Append(converged);
            sb.Append($" A total of {IterationsRequired} iterations were run and");
            sb.Append($" norm2(rhs - matrix * xSolution) / norm2(rhs - matrix * xInit) = {NormRatio}.");
            return sb.ToString();
        }
    }
}
