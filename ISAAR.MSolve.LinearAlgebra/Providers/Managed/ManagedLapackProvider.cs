using System;
using DotNumerics.LinearAlgebra.CSLapack;
using ISAAR.MSolve.LinearAlgebra.Commons;

//TODO: find a managed BLAS that supports the methods DotNumerics doesn't.
//TODO: In the custom LAPACK implementations, provide error checking for more than just the index where singularity, etc. is 
//      found. See LAPACK source for the checks. Error checking needs to be improved in general.
namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Provides managed C# implementations of the linear algebra operations defined by <see cref="ILapackProvider"/>. Uses the 
    /// library DotNumerics (see http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSLapack/Default.aspx) for the most
    /// part. For LAPACK subroutines not provided by DotNumerics, custom C# implementations are used instead. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ManagedLapackProvider : ILapackProvider
    {
        private static readonly DGELQF dgelqf = new DGELQF();
        private static readonly DGEQRF dgeqrf = new DGEQRF();
        private static readonly DORGLQ dorglq = new DORGLQ();
        private static readonly DORGQR dorgqr = new DORGQR();
        private static readonly DORMLQ dormlq = new DORMLQ();
        private static readonly DORMQR dormqr = new DORMQR();
        private static readonly DGETRF dgetrf = new DGETRF();
        private static readonly DGETRI dgetri = new DGETRI();
        private static readonly DGETRS dgetrs = new DGETRS();
        private static readonly DTRSM dtrsm = new DTRSM();

        internal static ManagedLapackProvider UniqueInstance { get; } = new ManagedLapackProvider();

        private ManagedLapackProvider() { } // private constructor for singleton pattern

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgelqf.aspx
        /// </summary>
        public void Dgelqf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dgelqf.Run(m, n, ref a, offsetA, ldA, ref tau, offsetTau, ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgeqrf.aspx
        /// </summary>
        public void Dgeqrf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => dgeqrf.Run(m, n, ref a, offsetA, ldA, ref tau, offsetTau, ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgetrf.aspx
        /// </summary>
        public void Dgetrf(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, ref int info)
            => dgetrf.Run(m, n, ref a, offsetA, ldA, ref ipiv, offsetIpiv, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgetri.aspx
        /// </summary>
        public void Dgetri(int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dgetri.Run(n, ref a, offsetA, ldA, ipiv, offsetIpiv, ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dgetrs.aspx
        /// </summary>
        public void Dgetrs(string transA, int n, int nRhs, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] b, int offsetB, int ldB, ref int info)
            => dgetrs.Run(transA, n, nRhs, a, offsetA, ldA, ipiv, offsetIpiv, ref b, offsetB, ldB, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dorglq.aspx
        /// </summary>
        public void Dorglq(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dorglq.Run(m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dorgqr.aspx
        /// </summary>
        public void Dorgqr(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => dorgqr.Run(m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dormlq.aspx
        /// </summary>
        public void Dormlq(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, 
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => dormlq.Run(side, transQ, m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref c, offsetC, ldC,
                ref work, offsetWork, lWork, ref info);

        /// <summary>
        /// See http://www.dotnumerics.com/NumericalLibraries/LinearAlgebra/CSharpCodeFiles/dormqr.aspx
        /// </summary>
        public void Dormqr(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau,
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => dormqr.Run(side, transQ, m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref c, offsetC, ldC,
                ref work, offsetWork, lWork, ref info);

        public void Dpotrf(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
        {
            try
            {
                if (IsUpper(uplo)) LapackImplementations.CholeskyUpperFullColMajor(n, a, offsetA, ldA, ref info);
                else LapackImplementations.CholeskyLowerFullColMajor(n, a, offsetA, ldA, ref info);
            }
            catch (ArgumentException)
            {
                info = -1;
            }
        }

        public void Dpotri(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
        {
            if (ldA != n) throw new NotImplementedException("dpotri only works for ldA=n");

            // Start with an identity matrix
            var inverse = new double[n * n];
            for (int i = 0; i < n; ++i) inverse[i * ldA + i] = 1.0;

            // Solve (L*L^T) * inverse = I or (U^T*U) * inverse = I
            int infoSolve = LapackUtilities.DefaultInfo;
            Dpotrs(uplo, n, n, a, offsetA, ldA, inverse, 0, n, ref infoSolve);

            // Copy the inverse matrix over the factorization
            Array.Copy(inverse, 0, a, offsetA, n * n);

            info = 0; //TODO add checks
        }

        public void Dpotrs(string uplo, int n, int nRhs, double[] a, int offsetA, int ldA, 
            double[] b, int offsetB, int ldB, ref int info)
        {
            try
            {
                if (IsUpper(uplo)) // A*X=B <=> U^T * (U*X) = B
                {
                    dtrsm.Run("L", "U", "T", "N", n, nRhs, 1.0, a, offsetA, ldA, ref b, offsetB, ldB); // B = U^T \ B
                    dtrsm.Run("L", "U", "N", "N", n, nRhs, 1.0, a, offsetA, ldA, ref b, offsetB, ldB); // B = U \ B
                }
                else // A*X=B <=> L * (L^T*X) = B
                {
                    dtrsm.Run("L", "L", "N", "N", n, nRhs, 1.0, a, offsetA, ldA, ref b, offsetB, ldB); // B = L \ B
                    dtrsm.Run("L", "L", "T", "N", n, nRhs, 1.0, a, offsetA, ldA, ref b, offsetB, ldB); // B = L^T \ B
                }
                info = 0; // TODO: needs more checks
            }
            catch (ArgumentException)
            {
                info = -1;
            }
        }

        public void Dpptrf(string uplo, int n, double[] a, int offsetA, ref int info)
        {
            try
            {
                if (IsUpper(uplo)) LapackImplementations.CholeskyUpperPackedColMajor(n, a, offsetA, ref info);
                else LapackImplementations.CholeskyLowerPackedColMajor(n, a, offsetA, ref info);
            }
            catch (ArgumentException)
            {
                info = -1;
            }
        }

        public void Dpptri(string uplo, int n, double[] a, int offsetA, ref int info)
        {
            try
            {
                // Start with an identity matrix
                var inverse = new double[n * n];
                for (int i = 0; i < n; ++i) inverse[i * n + i] = 1.0;

                // Solve (L*L^T) * inverse = I or (U^T*U) * inverse = I
                int infoSolve = LapackUtilities.DefaultInfo;
                Dpptrs(uplo, n, n, a, offsetA, inverse, 0, n, ref infoSolve);

                // Copy the inverse matrix over the factorization
                if (IsUpper(uplo)) Conversions.FullColMajorToPackedUpperColMajor(n, inverse, a, offsetA);
                else Conversions.FullColMajorToPackedLowerColMajor(n, inverse, a, offsetA);
                info = 0; // TODO: needs more checks
            }
            catch (ArgumentException)
            {
                info = -1;
            }
        }

        public void Dpptrs(string uplo, int n, int nRhs, double[] a, int offsetA, double[] b, int offsetB, int ldB, ref int info)
        {
            try
            {
                if (IsUpper(uplo))
                {
                    // Process each column separately
                    for (int i = 0; i < nRhs; ++i)
                    {
                        // b = U^T \ b
                        ManagedBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Upper, TransposeMatrix.Transpose,
                            DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * nRhs, 1);

                        // b = U \ b
                        ManagedBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Upper, TransposeMatrix.NoTranspose, 
                            DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * nRhs, 1);
                    }
                }
                else
                {
                    // Process each column separately
                    for (int i = 0; i < nRhs; ++i)
                    {
                        // b = L \ b
                        ManagedBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Lower, TransposeMatrix.NoTranspose, 
                            DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * n, 1);

                        // b = L^T \ b
                        ManagedBlasProvider.UniqueInstance.Dtpsv(StoredTriangle.Lower, TransposeMatrix.Transpose, 
                            DiagonalValues.NonUnit, n, a, offsetA, b, offsetB + i * n, 1);
                    }
                }
                info = 0; // TODO: needs more checks
            }
            catch (ArgumentException)
            {
                info = -1;
            }
        }

        private static bool IsUpper(string uplo)
        {
            if (uplo.Equals("L") || uplo.Equals("l")) return false;
            else if (uplo.Equals("U") || uplo.Equals("u")) return true;
            else throw new ArgumentException("Parameter uplo must be U, u, L or l");
        }
    }
}
