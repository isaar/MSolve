using System;
using System.Collections.Generic;
using System.Text;

//TODO: Code duplication between this and the CG version. They could be joined, by adopting the name ResDotPrecondResNew in both.
//TODO: Checking δnew <= ε^2 * δ0 as proposed by "Painless CG, section B3", is more efficient than normalizing and computing 
//      square roots. Also the order must be δnew <= ε (ε*δ0), since ε^2 could be very small and risk precision loss. 
//      However this idea is less flexible and does not fit well with logging sqrt(r*s)/sqrt(r0*s0) or whatever metric is used 
//      for converging. The performance gain is probably negligible, but it should be benchmarked, especially for problems with  
//      relatively fast matrix-vector and vector-vector operations.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Implements the most common way to check the convergence of CG: sqrt(r*inv(M)*r) / sqrt(r0*inv(M)*r0) &lt;= tolerance.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RegularPcgConvergence : IPcgResidualConvergence
    {
        private double denominator;

        /// <summary>
        /// See <see cref="IPcgResidualConvergence.EstimateResidualNormRatio(PcgAlgorithmBase)"/>
        /// </summary>
        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) => Math.Sqrt(pcg.ResDotPrecondRes) / denominator;

        /// <summary>
        /// See <see cref="IPcgResidualConvergence.Initialize(PcgAlgorithmBase)"/>
        /// </summary>
        public void Initialize(PcgAlgorithmBase pcg) => denominator = Math.Sqrt(pcg.ResDotPrecondRes); 
    }
}
