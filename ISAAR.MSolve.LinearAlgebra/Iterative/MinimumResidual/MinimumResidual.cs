using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: perhaps some vectors and operations are not needed or can be simplified 
//      in the version without preconditioning
//TODO: use invalid values as initial
//TODO: remove Matlab/Fortran notation from comments
//TODO: the orthogonalizer should not be private.
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
//TODO: many vectors are continuously allocated and then GCed. This can now be fixed.
namespace ISAAR.MSolve.LinearAlgebra.Iterative.MinimumResidual
{
    /// <summary>
    /// Implements the MINRES algorithm for solving an n-by-n system of linear equations: A*x = b, where A is symmetric and b  
    /// is a given vector of length n. A may be indefinite. It can also be singular, in which case the least squares problem is 
    /// solved instead. The MINRES method is presented by C. C. Paige, M. A. Saunders in 
    /// https://www.researchgate.net/publication/243578401_Solution_of_Sparse_Indefinite_Systems_of_Linear_Equations.
    /// Reorthogonalization is introduced by D. Maddix. For more information, including the Matlab scripts from which this code 
    /// is ported from, see http://web.stanford.edu/group/SOL/software/minres/.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class MinRes
    {
        private const double eps = double.Epsilon; // TODO: replace it in the code, when all else has been successfully ported

        private readonly bool checkMatrixSymmetricity;
        private readonly int maxIterations;
        private readonly int numStoredOrthogonalDirections;
        private readonly bool printIterations;
        private readonly double residualTolerance;

        /// <summary>
        /// Initializes a new instance of <see cref="MinRes"/> with the specified settings and convergence criteria.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">If norm2(b-A*x) / norm2(b-A*x0) &lt;= <paramref name="residualTolerance"/>, the 
        ///     algorithm will terminate, where x is the current solution vector and x0 the initial guess.</param>
        /// <param name="numStoredOrthogonalDirections">If &gt;0 local reorthogonalization will be used to improve convergence. 
        ///     However that requires extra memory equal to <paramref name="numStoredOrthogonalDirections"/> * matrixOrder. 
        ///     The author suggests the values 10, 20 for memory economy or 50, 100 if the memory requirements can be met.
        ///     It must not be greater than the order of the matrix though.</param>
        /// <param name="checkMatrixSymmetricity">If true, the matrix A of a linear system A*x=b will be checked to verify it is 
        ///     symmetric, which is safer but requires an additional matrix-vector multiplication.</param>
        /// <param name="printIterations">If true, the current solution vector x, the estimated condition number of the matrix A
        ///     and other statistics will be written to the console at each iteration of the algorithm, which will hinder 
        ///     performance.</param>
        public MinRes(int maxIterations, double residualTolerance, int numStoredOrthogonalDirections = 0,
             bool checkMatrixSymmetricity = false, bool printIterations = false)
        {
            this.maxIterations = maxIterations;
            this.residualTolerance = residualTolerance;
            this.numStoredOrthogonalDirections = numStoredOrthogonalDirections;
            this.checkMatrixSymmetricity = checkMatrixSymmetricity;
            this.printIterations = printIterations;
        }

        /// <summary>
        /// Solves the linear system (A - s*I) * x = b, where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/> 
        /// and s = <paramref name="shift"/>. If the matrix A - s*I is singular, it solves the same linear least squares problem.
        /// </summary>
        /// <param name="matrix">The matrix of the original linear system. It must be symmetric.</param>
        /// <param name="rhsVector">The right hand side vector of the original linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="shift">A scalar parameter that controls the deviation of (A - s*I) * x = b from A * x = b.</param>
        public (IVector solution, MinresStatistics stats) Solve(IMatrixView matrix, IVector rhsVector, double shift = 0.0)
            => SolveInternal(new ExplicitMatrixTransformation(matrix), rhsVector, null, shift);

        /// <summary>
        /// Solves the linear system (A - s*I) * x = b preconditioned by the matrix <paramref name="preconditioner"/>, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/> and s = <paramref name="shift"/>. 
        /// If the matrix A - s*I is singular, it solves the same linear least squares problem.
        /// </summary>
        /// <param name="matrix">The matrix of the original linear system. It must be symmetric.</param>
        /// <param name="rhsVector">The right hand side vector of the original linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="preconditioner">A preconditioner matrix that has the same dimensions as A, but must be symmetric 
        ///     positive definite, contrary to A.</param>
        /// <param name="shift">A scalar parameter that controls the deviation of (A - s*I) * x = b from A * x = b.</param>
        public (IVector solution, MinresStatistics stats) Solve(IMatrixView matrix, IVector rhsVector, 
            IPreconditioner preconditioner, double shift = 0.0)
            => SolveInternal(new ExplicitMatrixTransformation(matrix), rhsVector, preconditioner, shift);

        /// <summary>
        /// Solves the linear system (A - s*I) * x = b, where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/> 
        /// and s = <paramref name="shift"/>. If the matrix A - s*I is singular, it solves the same linear least squares problem.
        /// </summary>
        /// <param name="matrix">The matrix of the original linear system. It must be symmetric.</param>
        /// <param name="rhsVector">The right hand side vector of the original linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="shift">A scalar parameter that controls the deviation of (A - s*I) * x = b from A * x = b.</param>
        public (IVector solution, MinresStatistics stats) Solve(ILinearTransformation matrix, IVector rhsVector, 
            double shift = 0.0)
            => SolveInternal(matrix, rhsVector, null, shift);

        /// <summary>
        /// Solves the linear system (A - s*I) * x = b preconditioned by the matrix <paramref name="preconditioner"/>, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/> and s = <paramref name="shift"/>. 
        /// If the matrix A - s*I is singular, it solves the same linear least squares problem.
        /// </summary>
        /// <param name="matrix">The matrix of the original linear system. It must be symmetric.</param>
        /// <param name="rhsVector">The right hand side vector of the original linear system. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="preconditioner">A preconditioner matrix that has the same dimensions as A, but must be symmetric 
        ///     positive definite, contrary to A.</param>
        /// <param name="shift">A scalar parameter that controls the deviation of (A - s*I) * x = b from A * x = b.</param>
        public (IVector solution, MinresStatistics stats) Solve(ILinearTransformation matrix, IVector rhsVector, 
            IPreconditioner preconditioner,  double shift = 0.0) 
            => SolveInternal(matrix, rhsVector, preconditioner, shift);

        private (IVector solution, MinresStatistics stats) SolveInternal(ILinearTransformation A, IVector b, 
            IPreconditioner M, double shift)
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

            IVector y;
            if (M == null) y = b.Copy();
            else
            {
                y = b.CreateZeroVectorWithSameFormat();
                M.SolveLinearSystem(b, y);
            }
                
            IVector r1 = b.Copy(); // initial guess x = 0 initial residual
            double beta1 = b.DotProduct(y);

            /// If b = 0 exactly, stop with x = 0.            

            CheckDefinitePreconditioner(M, beta1);
            if (beta1 == 0.0)
            {
                var stats = new MinresStatistics { TerminationCause = 0, ResidualNorm = 0.0, IterationsRequired = 0 };
                return (x, stats);
            }
            beta1 = Math.Sqrt(beta1); //Normalize y to get v1 later.
            if (checkMatrixSymmetricity)
            {
                if (M == null) CheckSymmetricMatrix(A, y, r1);
                else
                {
                    CheckSymmetricPreconditioner(M, y, r1);

                    // WARNING: originally this had y instead of b. However, some precision is lost when multiplying inv(M)*b.
                    // Coupled with the absurdely small tolerance used by IsMatrixSymmetric(), the matrix is flagged as 
                    // non-symmetric even though it is. By using b instead of y, the check is identical to the one done in 
                    // the unpreconditioned MINRES version.
                    CheckSymmetricMatrix(A, b, r1);
                }
            }

            /// ------------------------------------------------- 
            /// Initialize other quantities.
            /// -------------------------------------------------

            double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0;
            double qrnorm = beta1, phibar = beta1, rhs1 = beta1, rhs2 = 0.0;
            double tnorm2 = 0.0, gmax = 0.0, gmin = double.MaxValue;
            double cs = -1.0, sn = 0.0;
            IVector w = Vector.CreateZero(n);
            IVector w2 = Vector.CreateZero(n);
            IVector r2 = r1.Copy();


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
                IVector v = y.Copy();
                v.ScaleIntoThis(s); // v = vk if P = I
                orthogonalizer.StoreDirection(v); // Store old v for local reorthogonalization of new v, if it is enabled. 

                ShiftedMatrixVectorMult(A, v, y, shift); // shift is 0 otherwise solving A - shift*I

                //WARNING: the following works if itn is initialized and updated as in the matlab script
                if (itn >= 2) y.AxpyIntoThis(r1, - beta / oldb); // normalization is the division r1 by oldb 

                double alfa = v.DotProduct(y);              // alphak
                y.AxpyIntoThis(r2, - alfa / beta); // normalization of r2/beta = v

                // v will be normalized through y later. 
                // This is explicit orthogonalizing it versus the previous localSize lanczos vectors.
                orthogonalizer.Reorthogonalize(y);

                r1 = r2; // r1 is unnormalized vold. No need to copy: r2 will point to another vector after this

                // r2 is unnormalized v. In the preconditioned version I must copy it, because in the next iteration v = y
                // and v.ScaleIntoThis() will corrupt r2.
                if (M == null) r2 = y.Copy();
                else
                {
                    r2 = y;
                    M.SolveLinearSystem(r2, y);
                }

                oldb = beta; // oldb = betak
                beta = r2.DotProduct(y);  // beta = betak+1^2
                CheckDefinitePreconditioner(M, beta);
                beta = Math.Sqrt(beta);
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
                IVector w1 = w2.Copy();
                w2 = w;
                // Do efficiently: w = (v - oldeps * w1 - delta * w2) * denom 
                w = v.Axpy(w1, - oldeps);
                w.AxpyIntoThis(w2 , - delta);
                w.ScaleIntoThis(denom);
                x.AxpyIntoThis(w, phi);

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

        private static void CheckDefinitePreconditioner(IPreconditioner M, double beta)
        {
            if (beta < 0)
            {
                Debug.Assert(M != null); // If not, the dot product beta1 must be >= 0
                throw new IndefiniteMatrixException("The preconditioner M is not positive definite.");
            }
        }

        private static void CheckSymmetricMatrix(ILinearTransformation A, IVector y, IVector r1)
        {
            IVector w = y.CreateZeroVectorWithSameFormat();
            A.Multiply(y, w);
            IVector r2 = y.CreateZeroVectorWithSameFormat();
            A.Multiply(w, r2);
            double s = w.DotProduct(w);
            double t = y.DotProduct(r2);
            double z = Math.Abs(s - t);
            double epsa = (s + eps) * Math.Pow(eps, 1.0/3.0);
            if (z > epsa) throw new AsymmetricMatrixException("The matrix or linear transformation A is not symmetric.");
        }

        private void CheckSymmetricPreconditioner(IPreconditioner M, IVector y, IVector r1)
        {
            IVector r2 = y.CreateZeroVectorWithSameFormat();
            M.SolveLinearSystem(y, r2);
            double s = y.DotProduct(y);
            double t = r1.DotProduct(r2);
            double z = Math.Abs(s - t);
            double epsa = (s + eps) * Math.Pow(eps, 1.0 / 3.0);
            if (z > epsa) throw new AsymmetricMatrixException("The preconditioner M is not symmetric.");
        }

        /// <summary>
        /// Calculates (A - shift * I) * v = A*v - shift*v
        /// </summary>
        private static void ShiftedMatrixVectorMult(ILinearTransformation matrix, IVectorView x, IVector y, double shift)
        {
            //TODO: this should just be implemented as a wrapping LinearTransformation
            if (shift == 0.0) matrix.Multiply(x, y);
            else
            {
                matrix.Multiply(x, y);
                y.AxpyIntoThis(x, -shift);
            }
        }

        private static void WriteIterationData(ILinearTransformation A, IVector b, double shift, IVector x, int iter,
            int maxIterations, double test1, double test2, double Anorm, double Acond, double gbar, double qrnorm, 
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
                    // Print true Arnorm. This works only if there is no preconditioning

                    // vv = b - (A - shift * I) * x
                    IVector vv = b.CreateZeroVectorWithSameFormat();
                    ShiftedMatrixVectorMult(A, x, vv, shift);
                    vv.LinearCombinationIntoThis(-1.0, b, 1.0);

                    // ww = (A - shift * I) * vv = "Ar"
                    IVector ww = b.CreateZeroVectorWithSameFormat();
                    ShiftedMatrixVectorMult(A, vv, ww, shift);
                    double trueArnorm = ww.Norm2();
                    Console.WriteLine();
                    Console.WriteLine($"Arnorm = {Anorm}  - True |Ar| = {trueArnorm}");
                }
            }
        }

        private class LocalReorthogonalizer
        {
            private readonly int localSize;
            private readonly LinkedList<IVector> store;

            internal LocalReorthogonalizer(int order, int localSize)
            {
                if ((localSize < 0) || (localSize > order)) throw new ArgumentException(
                    "The number of stored vectors must belong to [0, matrix.Order), but was " + localSize);
                this.localSize = localSize;
                this.store = new LinkedList<IVector>();
            }

            /// <summary>
            /// If orthogonalization is enabled, store a new direction <paramref name="v"/>.
            /// </summary>
            /// <param name="v"></param>
            internal void StoreDirection(IVector v)
            {
                if (localSize > 0) store.AddLast(v);
                if (store.Count > localSize) store.RemoveFirst();
            }

            /// <summary>
            /// Gram–Schmidt orthogonalization to the stored vectors, without normalization. 
            /// </summary>
            /// <param name="v">The vector to reorthogonalize. It will be overwritten with the result. If there are no stored 
            ///     vectors, it will remain unchanged without being copied.</param>
            internal void Reorthogonalize(IVector v)
            {
                // reorthogonalize 1 by 1
                foreach (var p in store)
                {
                    // we don't have to normalize since it is explicitly done in the code and so we don't need to redo it
                    v.AxpyIntoThis(p , - (v.DotProduct(p))); // orthogonalize to each stored vector
                }
            }
        }
    }
}
