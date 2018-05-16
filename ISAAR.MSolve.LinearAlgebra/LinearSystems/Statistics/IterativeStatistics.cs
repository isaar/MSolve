using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Add time measurements, flop measurements, etc. 
//TODO: Alternatively I could use a logger (pull observer) in the algorithm
//TODO: Needs a better name
//TODO: Each algorithm/author outputs something different. Once enough have been implemented/ported, find an appropriate design
//      to unify them. 
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics
{
    public class IterativeStatistics
    {
        public string AlgorithmName { get; set; }

        public bool HasConverged { get; set; }

        public int IterationsRequired { get; set; }

        /// <summary>
        /// lastResidual.Norm2() / firstResidualNorm2()
        /// </summary>
        public double NormRatio { get; set; }

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
