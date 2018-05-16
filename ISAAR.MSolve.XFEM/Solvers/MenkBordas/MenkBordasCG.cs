using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class MenkBordasCG
    {
        private readonly int maxIterations;
        private readonly double tolerance;

        public MenkBordasCG(int maxIterations, double tolerance)
        {
            this.maxIterations = maxIterations;
            this.tolerance = tolerance;
        }

        public (MenkBordasVector, IterativeStatistics) Solve(MenkBordasSystem sys)
        {
            sys.CheckDimensions();
            if ((sys.xs == null) || (sys.xe == null))
            {
                // 0) r = b - K * x
                MenkBordasVector r = sys.BuildRhsVector().Copy();
                var x = MenkBordasVector.CreateZeroWithSameDimensions(r);
                return SolveInternal(sys.BuildMatrix(), r, x);
            }
            else throw new NotImplementedException();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="x">The solution vector initialized appropriately. It will also be returned.</param>
        /// <returns></returns>
        private (MenkBordasVector x, IterativeStatistics stats) SolveInternal(MenkBordasMatrix A, MenkBordasVector r,
            MenkBordasVector x)
        {
            // 1) p[0] = r[0]
            MenkBordasVector p = r.Copy();

            // 2) dot[0] = r[0]^T * r[0]
            double rDot = r.DotProduct(r);

            // 3) norm2[0] = sqrt(dot(r[0]))
            var convergence = new Convergence(tolerance, r);

            IterativeStatistics statistics;
            for (int t = 0; t < maxIterations; ++t)
            {
                // 4) Ap = A * p
                MenkBordasVector Ap = A.MultiplyRight(p);

                // 5) a = dot[t] / p^T*K*p
                double a = rDot / p.DotProduct(Ap);

                // 6) x = x + a * p
                x.AxpyIntoThis(a, p);

                // 7) r = r - a * Ap 
                r.AxpyIntoThis(-a, Ap);

                // 8) dot[t+1] = r[t+1]^T * r[t+1] 
                double rDotNext = r.DotProduct(r);

                // 9) convergence check: norm2[t] / norm2[0] < tol 
                if (convergence.HasConvergedGlobally(r))
                {
                    statistics = new IterativeStatistics
                    {
                        AlgorithmName = "Conjugate Gradient for Menk-Bordas system",
                        HasConverged = true,
                        IterationsRequired = t + 1,
                        NormRatio = convergence.CurrentResNormRatio(rDot)
                    };
                    return (x, statistics);
                }

                // 10) b = dot(r[t+1]) / dot(r[t])
                double b = rDotNext / rDot;

                // 11) p = r + p * b
                p = r.Axpy(b, p);

                // 12) Swap the dot product to avoid recomputing it.
                rDot = rDotNext;
            }

            statistics = new IterativeStatistics
            {
                AlgorithmName = "Conjugate Gradient",
                HasConverged = false,
                IterationsRequired = maxIterations,
                NormRatio = convergence.CurrentResNormRatio(r.DotProduct(r)) // Actually I think I need the previous r here.
            };
            return (x, statistics);
        }

        private class Convergence
        {
            private readonly double tolerance;
            private readonly double[] rDotsInit;
            private readonly double[] rNorms2Init;
            private readonly double totalNorm2Init;

            public Convergence(double tolerance, MenkBordasVector rInit)
            {
                this.tolerance = tolerance;
                rDotsInit = rInit.IndividualDots(rInit);
                rNorms2Init = new double[rDotsInit.Length];
                totalNorm2Init = 0;
                for (int i = 0; i < rDotsInit.Length; ++i)
                {
                    rNorms2Init[i] = Math.Sqrt(rDotsInit[i]);
                    totalNorm2Init += rDotsInit[i];
                }
                totalNorm2Init = Math.Sqrt(totalNorm2Init);
            }

            public bool HasConvergedIndividually(MenkBordasVector r)
            {
                double[] rNorms2 = r.IndividualNorms2();
                for (int i = 0; i < rNorms2.Length; ++i)
                {
                    if (rNorms2[i] / rNorms2Init[i] > tolerance) return false;
                }
                return true;
            }

            public bool HasConvergedGlobally(MenkBordasVector r)
            {
                return Math.Sqrt(r.DotProduct(r)) / totalNorm2Init <= tolerance;
            }

            public double CurrentResNormRatio(double rDot)
            {
                return Math.Sqrt(rDot) / totalNorm2Init;
            }
        }
    }
}
