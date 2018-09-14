using System;
using System.Text;

//TODO: not too thrilled by this style of reporting. Design it in an OOP way with excpetions, etc.
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.MinRes
{
    /// <summary>
    /// Data Transfer Object that collects information about the execution and convergence when solving a linear system with 
    /// the Minimum Redidual algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MinresStatistics
    {
        /// <summary>
        /// The number of iterations that were run before the algorithm terminated. 
        /// </summary>
        public int IterationsRequired { get; set; } = int.MinValue;

        /// <summary>
        /// A code that describes why the algorithm has terminated. For a more readable form, use
        /// <see cref="TerminationMessage"/> instead.
        /// </summary>
        public int TerminationCause { get; set; } = int.MinValue; // TODO: make an enum for this, that will also be used in the algorithm

        /// <summary>
        /// A message that describes why the algorithm has terminated.
        /// </summary>
        public string TerminationMessage
        {
            get
            {
                switch (TerminationCause)
                {
                    case -1:
                        return "beta2 = 0. If M = I, b and x are eigenvectors.";
                    case 0:
                        return "beta1 = 0. The exact solution is x = 0.";
                    case 1:
                        return "A solution to Ax = b was found, given the tolerance in norm2(residual).";
                    case 2:
                        return "A least-squares solution was found, given the tolerance in norm2(residual).";
                    case 3:
                        return "Reasonable accuracy achieved, given eps=double.Epsilon";
                    case 4:
                        return "x has converged to an eigenvector.";
                    //case 5: //TODO: find out why this is never ouput. Instead when Acond > HUGE NUMBER case 4 is returned! 
                    //    return "Acond has exceeded 0.1/double.Epsilon."; //TODO: This means divergence right? Should I throw an exception?
                    case 6:
                        return "The iteration limit was reached.";
                    case int.MinValue:
                        throw new ArgumentNullException("The termination cause has not been set.");
                    default:
                        throw new ArgumentException("Invalid termination cause.");
                }
            }
        }

        /// <summary>
        /// Estimate for cond(A). In the current version of MINRES, it is computed from the diagonals of R in the factorization
        ///  of the lower Hessenberg matrix, Q* H = R, where H is the tridiagonal matrix from Lanczos with one extra row, 
        /// beta(k + 1) e_k^T.
        /// </summary>
        public double MatrixCondition { get; set; } = double.MinValue;

        /// <summary>
        /// Estimate for the induced norm2 of the matrix A.
        /// </summary>
        public double MatrixNorm { get; set; } = double.MinValue;

        /// <summary>
        /// Estimate for norm2(A * (b-A*x) ) at the previous to last iteration. Note that it lags one iteration 
        /// behind <see cref="ResidualNorm"/>.
        /// </summary>
        public double MatrixTimesResidualNorm { get; set; } = double.MinValue;

        /// <summary>
        /// Estimate for norm2(b-A*x) at the last iteration.
        /// </summary>
        public double ResidualNorm { get; set; } = double.MinValue;

        /// <summary>
        /// Value of norm2(x) at the last iteration.
        /// </summary>
        public double YNorm { get; set; } = double.MinValue;

        /// <summary>
        /// Reports the accumulated data of this <see cref="MinresStatistics"/> instance. 
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(TerminationMessage);
            sb.AppendLine($"Iterations required = {IterationsRequired}");
            sb.AppendLine($"Estimated matrix norm = {MatrixNorm}");
            sb.AppendLine($"Estimated matrix condition = {MatrixCondition}");
            sb.AppendLine($"Norm2 of last residual vector = {ResidualNorm}");
            sb.AppendLine($"Norm2 of last 'y' vector = {YNorm}");
            sb.AppendLine($"Norm2 of previous than last matrix * residual vector = { MatrixTimesResidualNorm}");
            
            return sb.ToString();
        }
    }
}
