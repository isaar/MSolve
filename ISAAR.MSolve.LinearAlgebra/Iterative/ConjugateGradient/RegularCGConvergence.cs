using System;
using System.Collections.Generic;
using System.Text;

//TODO: Checking δnew <= ε^2 * δ0 as proposed by "Painless CG, section B2", is more efficient than normalizing and computing 
//      square roots. Also the order must be δnew <= ε (ε*δ0), since ε^2 could be very small and risk precision loss. 
//      However this idea is less flexible and does not fit well with logging |r|/|r0| or whatever metric is used for converging. 
//      The performance gain is probably negligible, but it should be benchmarked, especially for problems with relatively fast 
//      matrix-vector and vector-vector operations.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the most common way to check the convergence of CG: norm2(r) / norm2(r0) &lt;= tolerance.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RegularCGConvergence : ICGResidualConvergence
    {
        private double residualNormInitial;

        /// <summary>
        /// See <see cref="ICGResidualConvergence.EstimateResidualNormRatio(CGAlgorithm)"/>
        /// </summary>
        public double EstimateResidualNormRatio(CGAlgorithm cg) => Math.Sqrt(cg.ResDotRes) / residualNormInitial;

        /// <summary>
        /// See <see cref="ICGResidualConvergence.Initialize(CGAlgorithm)"/>
        /// </summary>
        public void Initialize(CGAlgorithm cg) => this.residualNormInitial = Math.Sqrt(cg.ResDotRes); 
    }
}
