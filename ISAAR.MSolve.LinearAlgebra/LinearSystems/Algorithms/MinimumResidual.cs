using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: perhaps some vectors and operations are not needed or can be simplified 
//      in the version without preconditioning
//TODO: use invalid values as initial
//TODO: remove Matlab/Fortran notation from comments
//TODO: the orthogonalizer should not be private.
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms
{
    /// <summary>
    /// MINRES algorithm for solving an n-by-n system of linear equations: A*x = b, where A is symmetric and b is a given vector 
    /// of length n. A may be indefinite. It can also be singular, in which case the least squares problem is solved instead.
    /// The MINRES method is presented by C. C. Paige, M. A. Saunders in 
    /// https://www.researchgate.net/publication/243578401_Solution_of_Sparse_Indefinite_Systems_of_Linear_Equations.
    /// Reorthogonalization is introduced by D. Maddix. For more information, including the Matlab scripts from which this code 
    /// is ported from, see http://web.stanford.edu/group/SOL/software/minres/.
    /// </summary>
    public class MinimumResidual
    {
        private const double eps = double.Epsilon; // TODO: replace it in the code, when all else has been successfully ported

        private readonly bool checkMatrixSymmetricity;
        private readonly int maxIterations;
        private readonly int numStoredOrthogonalDirections;
        private readonly bool printIterations;
        private readonly double residualTolerance;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="maxIterations"></param>
        /// <param name="residualTolerance"></param>
        /// <param name="numStoredOrthogonalDirections">If &gt;0 local reorthogonalization will be used to improve convergence. 
        ///     However that requires extra memory equal to <paramref name="numStoredOrthogonalDirections"/> * matrixOrder. 
        ///     The author suggests the values 10, 20 for memory economy or 50, 100 if the memory requirements can be met.
        ///     It must not be greater than the order of the matrix though.</param>
        /// <param name="checkMatrixSymmetricity"></param>
        /// <param name="printIterations"></param>
        public MinimumResidual(int maxIterations, double residualTolerance, int numStoredOrthogonalDirections = 0, 
            bool checkMatrixSymmetricity = false, bool printIterations = false)
        {
            this.maxIterations = maxIterations;
            this.residualTolerance = residualTolerance;
            this.numStoredOrthogonalDirections = numStoredOrthogonalDirections;
            this.checkMatrixSymmetricity = checkMatrixSymmetricity;
            this.printIterations = printIterations;
        }

        public (Vector, MinresStatistics) Solve(IMatrixView A, Vector b, double shift = 0.0)
        {
            /// Initialize.

            int n = b.Length;
            int istop = 0; // TODO: not sure if needed
            int itn = 0;
            double Anorm = 0.0;
            double Acond = 0.0;
            double rnorm = 0.0;
            double ynorm = 0.0;
            var x = Vector.CreateZero(n);
            var orthogonalizer = new LocalReorthogonalizer(n, numStoredOrthogonalDirections);

            /// -------------------------------------------------
            /// Set up y and v for the first Lanczos vector v1.
            /// y = beta1 P' v1,  where  P = C**(-1).
            /// v is really P' v1.
            /// -------------------------------------------------

            Vector y = b.Copy();
            Vector r1 = b.Copy(); // initial guess x = 0 initial residual
            double beta1 = b * y;

            /// If b = 0 exactly, stop with x = 0.            

            if (beta1 == 0.0)
            {
                var stats = new MinresStatistics { TerminationCause = 0, ResidualNorm = 0.0, IterationsRequired = 0 };
                return (x, stats);
            }
            beta1 = Math.Sqrt(beta1); //Normalize y to get v1 later.
            if (checkMatrixSymmetricity) CheckSymmetricMatrix(A, y, r1);

            /// ------------------------------------------------- 
            /// Initialize other quantities.
            /// -------------------------------------------------

            double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0;
            double qrnorm = beta1, phibar = beta1, rhs1 = beta1, rhs2 = 0.0;
            double tnorm2 = 0.0, gmax = 0.0, gmin = double.MaxValue;
            double cs = -1.0, sn = 0.0;
            var w = Vector.CreateZero(n);
            var w2 = Vector.CreateZero(n);
            Vector r2 = r1.Copy();


            /// ------------------------------------------------- 
            /// Main iteration loop.
            /// -------------------------------------------------
            while (itn < maxIterations) // k = itn = 1 first time through
            {
                ++itn;

                /// -----------------------------------------------------------------
                /// Obtain quantities for the next Lanczos vector vk + 1, k = 1, 2,...
                ///The general iteration is similar to the case k = 1 with v0 = 0:
                ///
                /// p1 = Operator * v1 - beta1 * v0,
                /// alpha1 = v1'p1,
                /// q2 = p2 - alpha1 * v1,   // I think this comment was corrected to q2 = p1 - alpha1 * v1 in a later script
                /// beta2 ^ 2 = q2'q2,
                /// v2 = (1 / beta2) q2.
                ///
                /// Again, y = betak P vk, where  P = C * *(-1).
                /// ....more description needed.
                /// -----------------------------------------------------------------

                double s = 1.0 / beta;    // Normalize previous vector(in y).
                Vector v = y; // No need to copy: y will point to another vector shortly after this
                v.ScaleIntoThis(s); // v = vk if P = I
                orthogonalizer.StoreDirection(v); // Store old v for local reorthogonalization of new v, if it is enabled. 

                y = ShiftedMatrixVectorMult(A, v, shift); // shift is 0 otherwise solving A - shift*I

                //WARNING: the following works if itn is initialized and updated as in the matlab script
                if (itn >= 2) y.AxpyIntoThis(-beta / oldb, r1); // normalization is the division r1 by oldb 

                double alfa = v * y;              // alphak
                y.AxpyIntoThis(-alfa / beta, r2); // normalization of r2/beta = v

                // v will be normalized through y later. 
                // This is explicit orthogonalizing it versus the previous localSize lanczos vectors.
                orthogonalizer.Reorthogonalize(y);


                r1 = r2; // r1 is unnormalized vold. No need to copy: r2 will point to another vector after this
                r2 = y.Copy(); // r2 is unnormalized v. In the preconditioned version I must copy it, because in the next iteration v = y and v.ScaleIntoThis() will corrupt r2.
                oldb = beta; // oldb = betak
                beta = Math.Sqrt(r2 * y);  // beta = betak+1^2
                tnorm2 += alfa * alfa + oldb * oldb + beta * beta;

                if (itn == 1) //Initialize a few things.
                {
                    if (beta / beta1 <= 10.0 * eps) //beta2 = 0 or ~0
                    {
                        istop = -1; // Terminate later.
                    }
                }

                /// Apply previous rotation Qk-1 to get
                ///   [deltak epslnk + 1] = [cs  sn][dbark    0]
                ///   [gbar k dbar k + 1]   [sn - cs][alfak betak + 1].

                double oldeps = epsln;
                double delta = cs * dbar + sn * alfa; // delta1 = 0         deltak
                double gbar = sn * dbar - cs * alfa;  // gbar 1 = alfa1     gbar k
                double gbarSquared = gbar * gbar;
                epsln = sn * beta;                    // epsln2 = 0         epslnk + 1
                dbar = -cs * beta;                    // dbar 2 = beta2     dbar k+1
                double root = Math.Sqrt(gbarSquared + dbar * dbar);
                double Arnorm = phibar * root;        // || Ar{ k - 1}||

                /// Compute the next plane rotation Qk

                double gamma = Math.Sqrt(gbarSquared + beta * beta); // gammak
                if (gamma < eps) gamma = eps;
                cs = gbar / gamma;                                   // ck
                sn = beta / gamma;                                   // sk
                double phi = cs * phibar;                                   // phik
                phibar = sn * phibar;                                // phibark + 1

                /// Update  x.

                double denom = 1.0 / gamma;
                Vector w1 = w2.Copy();
                w2 = w;
                // Do efficiently: w = (v - oldeps * w1 - delta * w2) * denom 
                w = v.Axpy(-oldeps, w1);
                w.AxpyIntoThis(-delta, w2);
                w.ScaleIntoThis(denom);
                x.AxpyIntoThis(phi, w);

                /// Go round again.

                if (gmax < gamma) gmax = gamma;
                if (gmin > gamma) gmin = gamma;
                double z = rhs1 / gamma;
                rhs1 = rhs2 - delta * z;
                rhs2 = -epsln * z;

                /// Estimate various norms.

                Anorm = Math.Sqrt(tnorm2);
                ynorm = x.Norm2();
                double epsa = Anorm * eps;
                double epsx = Anorm * ynorm * eps;
                double epsr = Anorm * ynorm * residualTolerance;
                double diag = gbar;
                if (diag == 0) diag = epsa;

                qrnorm = phibar;
                rnorm = qrnorm;
                double test1 = rnorm / (Anorm * ynorm);    //  || r || / (|| A || || x ||)
                double test2 = root / Anorm;               // || Ar{ k - 1}|| / (|| A || || r_{ k - 1}||)


                /// Estimate  cond(A).
                /// In this version we look at the diagonals of R  in the
                /// factorization of the lower Hessenberg matrix,  Q* H = R,
                /// where H is the tridiagonal matrix from Lanczos with one
                /// extra row, beta(k + 1) e_k ^ T.

                Acond = gmax / gmin;

                /// See if any of the stopping criteria are satisfied.
                /// In rare cases, istop is already - 1 from above (Abar = const* I).

                if (istop == 0)
                {
                    double t1 = 1.0 + test1;      // These tests work if residualTolerance < double.Epsilon. 
                    double t2 = 1.0 + test2;      //TODO: then test that directly, like the last ones. Geez
                    if (t2 <= 1.0) istop = 2;
                    if (t1 <= 1.0) istop = 1;

                    if (itn >= maxIterations) istop = 6; 
                    if (Acond >= 0.1 / eps) istop = 4; // TODO: shouldn't this be case istop=5?
                    if (epsx >= beta1) istop = 3;
                    //%if rnorm <= epsx   , istop = 2; end //They were commented out in the original code
                    //%if rnorm <= epsr   , istop = 1; end //They were commented out in the original code
                    if (test2 <= residualTolerance) istop = 2;
                    if (test1 <= residualTolerance) istop = 1;
                }

                if (printIterations) WriteIterationData(A, b, shift, x, itn, maxIterations, test1, test2, Anorm, Acond, 
                    gbar, qrnorm, istop, epsx, epsr);

                if (istop != 0)
                {
                    var stats = new MinresStatistics
                    {
                        IterationsRequired = itn, //WARNING: take care if you change itn to be a for loop variable or sth else
                        TerminationCause = istop,
                        MatrixCondition = Acond,
                        MatrixNorm = Anorm,
                        MatrixTimesResidualNorm = Arnorm,
                        ResidualNorm = rnorm,
                        YNorm = ynorm
                    };
                    return (x, stats);
                }
            }
            throw new Exception("Should not have reached here");
        }

        internal static void CheckSymmetricMatrix(IMatrixView A, Vector y, Vector r1)
        {
            Vector w = A.MultiplyRight(y, false);
            Vector r2 = A.MultiplyRight(w, false);
            double s = w * w;
            double t = y * r2;
            double z = Math.Abs(s - t);
            double epsa = (s + eps) * Math.Pow(eps, 1.0/3.0);
            if (z > epsa) throw new AsymmetricMatrixException("The matrix or linear transformation A is not symmetric.");
        }

        /// <summary>
        /// Calculates (A - shift * v) = A*v - shift*v
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="vector"></param>
        /// <param name="shift"></param>
        /// <returns></returns>
        internal static Vector ShiftedMatrixVectorMult(IMatrixView matrix, Vector vector, double shift)
        {
            if (shift == 0.0) return matrix.MultiplyRight(vector, false);
            else
            {
                Vector result = matrix.MultiplyRight(vector, false);
                result.AxpyIntoThis(shift, vector);
                return result;
            }
        }

        internal static void WriteIterationData(IMatrixView A, Vector b, double shift, Vector x, int iter, int maxIterations,
            double test1, double test2, double Anorm, double Acond, double gbar, double qrnorm, 
            int istop, double epsx, double epsr)
        {
            bool printThisOnce = false;
            if (x.Length < 40) printThisOnce = true;
            if (iter < 10) printThisOnce = true;
            if (iter >= maxIterations - 11) printThisOnce = true;
            if (iter / 10 == 0) printThisOnce = true;
            if (qrnorm <= 10 * epsx) printThisOnce = true;
            if (qrnorm <= 10 * epsr) printThisOnce = true;
            if (Acond <= 1e-2 / eps) printThisOnce = true;
            if (istop != 0) printThisOnce = true;

            if (printThisOnce)
            {
                if (iter % 10 == 0) Console.WriteLine();
                Console.WriteLine($"Iteration {iter}: x[0] = {x[0]}, Compatible = {test1}, LS = {test2},"
                    + $" |A| = {Anorm}, cond(A) = {Acond}, gbar/|A| = {gbar / Anorm}");

                bool debug = false;
                if (debug)
                {
                    // Print true Arnorm. This works only if no preconditioning
                    Vector vv = b - ShiftedMatrixVectorMult(A, x, shift);   // vv = b - (A - shift * I) * x
                    Vector ww = ShiftedMatrixVectorMult(A, vv, shift);      // ww = (A - shift * I) * vv = "Ar"
                    double trueArnorm = ww.Norm2();
                    Console.WriteLine();
                    Console.WriteLine($"Arnorm = {Anorm}  - True |Ar| = {trueArnorm}");
                }
            }
        }

        class LocalReorthogonalizer
        {
            private readonly int localSize;
            private readonly LinkedList<Vector> store;

            internal LocalReorthogonalizer(int order, int localSize)
            {
                if ((localSize < 0) || (localSize > order)) throw new ArgumentException(
                    "The number of stored vectors must belong to [0, matrix.Order), but was " + localSize);
                this.localSize = localSize;
                this.store = new LinkedList<Vector>();
            }

            /// <summary>
            /// If orthogonalization is enabled, store a new direction <paramref name="v"/>.
            /// </summary>
            /// <param name="v"></param>
            internal void StoreDirection(Vector v)
            {
                if (localSize > 0) store.AddLast(v);
                if (store.Count > localSize) store.RemoveFirst();
            }

            /// <summary>
            /// Gram–Schmidt orthogonalization to the stored vectors, without normalization. 
            /// </summary>
            /// <param name="v">The vector to reorthogonalize. It will be overwritten with the result. If there are no stored 
            ///     vectors, it will remain unchanged without being copied.</param>
            internal void Reorthogonalize(Vector v)
            {
                // reorthogonalize 1 by 1
                foreach (var p in store)
                {
                    // we don't have to normalize since it is explicitly done in the code and so we don't need to redo it
                    v.AxpyIntoThis(-(v * p), p); // orthogonalize to each stored vector
                }
            }
        }
    }
}
