using System;
using static ISAAR.MSolve.LinearAlgebra.Providers.LapackUtilities;

//TODO: provide versions of these methods where the user can provide the work arrays.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Simplifies the use of LAPACK (see <see cref="ILapackProvider"/>) linear algebra operations that concern the solution of
    /// systems of linear equations with double precision arithmetic. Such simplifications are error checking, 
    /// handling workspace arrays, enums instead of string arguments etc. This facade is meant to provide a managed 
    /// object-oriented alternative the LAPACKE library used in C.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class LapackLinearEquationsFacade
    {
        internal LapackLinearEquationsFacade(ILapackProvider provider)
        {
            this.Provider = provider;
        }

        internal ILapackProvider Provider { get; }

        /// <summary>
        /// If the return value is a non-negative int i, then the pivot entry U[i, i] is 0. The factor U is singular, and 
        /// division by zero will occur if it is used to solve a system of equations with 
        /// <see cref="Dgetrs(TransposeMatrix, int, int, double[], int, int, int[], int, double[], int, int)"/>. Matrix inversion
        /// using <see cref="Dgetri(int, double[], int, int, int[], int)"/> will also fail.
        /// </summary>
        internal int Dgetrf(int numRowsA, int numColsA, double[] matrixA, int offsetA, int leadingDimA, 
            int[] rowExchangesP, int offsetP, double pivotTolerance)
        { 
            int info = DefaultInfo;
            Provider.Dgetrf(numRowsA, numColsA, matrixA, offsetA, leadingDimA, rowExchangesP, offsetP, ref info);

            if (info < 0)
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -4) info = -5;
                else if (info == -5) info = -6;
                ProcessNegativeInfo(info);
            }

            return FindZeroPivotFromInfoLU(numColsA, leadingDimA, matrixA, offsetA, pivotTolerance, info);
        }

        /// <summary>
        /// If the return value is a non-negative int i, then the pivot entry U[i, i] is 0. The matrix A is singular, and its
        /// inverse could not be computed.
        /// </summary>
        internal int Dgetri(int orderA, double[] factorizedMatrixA, int offsetA, int leadingDimA, 
            int[] rowExchangesP, int offsetP, double pivotTolerance)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dgetri(
                orderA, factorizedMatrixA, offsetA, leadingDimA, rowExchangesP, offsetP, work, offsetWork, lWork, ref info));

            if (info < 0)
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -3) info = -3;
                else if (info == -4) info = -5;
                ProcessNegativeInfo(info);
            }

            return FindZeroPivotFromInfoLU(orderA, leadingDimA, factorizedMatrixA, offsetA, pivotTolerance, info);
        }

        //internal void DgetriWithoutCheckingPivot(int orderA, double[] factorizedMatrixA, int offsetA, int leadingDimA,
        //    int[] rowExchangesP, int offsetP)
        //{
        //    int info = DefaultInfo;
        //    QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dgetri(
        //        orderA, factorizedMatrixA, offsetA, leadingDimA, rowExchangesP, offsetP, work, offsetWork, lWork, ref info));
        //    CheckNegativeInfo(info);
        //}

        internal void Dgetrs(TransposeMatrix transposeA, int orderA, int numRhs, double[] factorizedMatrixA, int offsetA,
            int leadingDimA, int[] rowExchangesP, int offsetP, double[] rhsB, int offsetB, int leadingDimB)
        {
            int info = DefaultInfo;
            Provider.Dgetrs(transposeA.Translate(), orderA, numRhs, factorizedMatrixA, offsetA, leadingDimA, 
                rowExchangesP, offsetP, rhsB, offsetB, leadingDimB, ref info);

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -5) info = -6;
                else if (info == -6) info = -7;
                else if (info == -7) info = -9;
                else if (info == -8) info = -11;
                ProcessNegativeInfo(info);
            }
        }

        /// <summary>
        /// If the return value is a non-negative int i, then the leading minor of order i (and the matrix A) is not positive
        /// definite and the factorization could not be completed.
        /// </summary>
        internal int Dpotrf(StoredTriangle triangle, int orderA, double[] matrixA, int offsetA, int leadingDimA)
        {
            int info = DefaultInfo;
            Provider.Dpotrf(triangle.Translate(), orderA, matrixA, offsetA, leadingDimA, ref info);

            if (info < 0)
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -4) info = -5;
                ProcessNegativeInfo(info);
            }

            //TODO: is there the same problem as LU, i.e. can the last minor not be positive definite, while the info = 0? 
            if (info > 0) return info - 1;
            else return int.MinValue;
        }

        /// <summary>
        /// If the return value is a non-negative int i, then the entry (i, i) of the factor U or L is 0 and the inverse could 
        /// not be computed.
        /// </summary>
        internal int Dpotri(StoredTriangle triangle, int orderA, double[] factorizedMatrixA, int offsetA, int leadingDimA)
        {
            int info = DefaultInfo;
            Provider.Dpotri(triangle.Translate(), orderA, factorizedMatrixA, offsetA, leadingDimA, ref info);

            if (info < 0)
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -4) info = -5;
                ProcessNegativeInfo(info);
            }

            //TODO: is there the same problem as LU, i.e. can the last minor not be positive definite, while the info = 0? 
            if (info > 0) return info - 1;
            else return int.MinValue;
        }

        internal void Dpotrs(StoredTriangle triangle, int orderA, int numRhs, double[] factorizedMatrixA, int offsetA,
            int leadingDimA, double[] rhsB, int offsetB, int leadingDimB)
        {
            int info = DefaultInfo;
            Provider.Dpotrs(triangle.Translate(), orderA, numRhs, factorizedMatrixA, offsetA, leadingDimA, 
                rhsB, offsetB, leadingDimB, ref info);

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -5) info = -6;
                else if (info == -6) info = -7;
                else if (info == -7) info = -9;
                ProcessNegativeInfo(info);
            }
        }

        /// <summary>
        /// If the return value is a non-negative int i, then the leading minor of order i (and the matrix A) is not positive
        /// definite and the factorization could not been completed.
        /// </summary>
        internal int Dpptrf(StoredTriangle triangle, int orderA, double[] matrixA, int offsetA)
        {
            int info = DefaultInfo;
            Provider.Dpptrf(triangle.Translate(), orderA, matrixA, offsetA, ref info);

            if (info < 0) ProcessNegativeInfo(info);

            //TODO: is there the same problem as LU, i.e. can the last minor not be positive definite, while the info = 0? 
            if (info > 0) return info - 1;
            else return int.MinValue;
        }

        /// <summary>
        /// If the return value is a non-negative int i, then the entry (i, i) of the factor U or L is 0 and the inverse could 
        /// not be computed.
        /// </summary>
        internal int Dpptri(StoredTriangle triangle, int orderA, double[] factorizedMatrixA, int offsetA)
        {
            int info = DefaultInfo;
            Provider.Dpptri(triangle.Translate(), orderA, factorizedMatrixA, offsetA, ref info);

            if (info < 0) ProcessNegativeInfo(info);

            //TODO: is there the same problem as LU, i.e. can the last minor not be positive definite, while the info = 0? 
            if (info > 0) return info - 1;
            else return int.MinValue;
        }

        internal void Dpptrs(StoredTriangle triangle, int orderA, int numRhs, double[] factorizedMatrixA, int offsetA, 
            double[] rhsB, int offsetB, int leadingDimB)
        {
            int info = DefaultInfo;
            Provider.Dpptrs(triangle.Translate(), orderA, numRhs, factorizedMatrixA, offsetA, 
                rhsB, offsetB, leadingDimB, ref info);

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -5) info = -6;
                else if (info == -6) info = -8;
                ProcessNegativeInfo(info);
            }
        }

        //private static int FindIndefiniteMinorFromInfoCholesky(int infoCholesky)
        //{
        //    //TODO: is there the same problem as LU, i.e. can the last minor not be positive definite, while the info = 0? 
        //    if (infoCholesky > 0) return infoCholesky - 1;
        //    else return int.MinValue;
        //}

        private static int FindZeroPivotFromInfoLU(int numColsA, int leadingDimA, double[] factorizedMatrixA, int offsetA,
            double pivotTolerance, int infoLU)
        {
            int firstZeroPivot = int.MinValue;
            if (infoLU == 0) // Supposedly everything went ok
            {
                //TODO: not sure about the use of numColsA in the next.
                if (Math.Abs(factorizedMatrixA[offsetA + numColsA * leadingDimA - 1]) <= pivotTolerance)
                {
                    // False Negative: info = 0, but LAPACK doesn't check the last diagonal entry!
                    firstZeroPivot =  numColsA - 1;
                }
            }
            else if (infoLU > 0) firstZeroPivot = infoLU - 1;
            return firstZeroPivot;
        }
    }
}
