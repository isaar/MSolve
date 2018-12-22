using static ISAAR.MSolve.LinearAlgebra.Providers.LapackUtilities;

//TODO: provide versions of these methods where the user can provide the work arrays.
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Simplifies the use of LAPACK (see <see cref="ILapackProvider"/>) linear algebra operations that concern the solution of
    /// least squares and minimum norm problems with double precision arithmetic. Such simplifications are error checking, 
    /// handling workspace arrays, enums instead of string arguments etc. This facade is meant to provide a managed 
    /// object-oriented alternative the LAPACKE library used in C.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class LapackLeastSquaresFacadeDouble
    {
        internal LapackLeastSquaresFacadeDouble(ILapackProvider provider)
        {
            this.Provider = provider;
        }

        internal ILapackProvider Provider { get; }

        internal void Dgelqf(int numRowsA, int numColsA, double[] matrixA, int offsetA, int leadingDimA,
            double[] reflectorScalarsT, int offsetT)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dgelqf(
                numRowsA, numColsA, matrixA, offsetA, leadingDimA, reflectorScalarsT, offsetT,
                work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -4) info = -5;
                else if (info == -5) info = -6;
                ProcessNegativeInfo(info);
            }
        }

        internal void Dgeqrf(int numRowsA, int numColsA, double[] matrixA, int offsetA, int leadingDimA,
            double[] reflectorScalarsT, int offsetT)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dgeqrf(
                numRowsA, numColsA, matrixA, offsetA, leadingDimA, reflectorScalarsT, offsetT,
                work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -4) info = -5;
                else if (info == -5) info = -6;
                ProcessNegativeInfo(info);
            }
        }

        internal void Dorglq(int numRowsQ, int numColsQ, int numReflectors, double[] matrixQ, int offsetQ, int leadingDimQ,
            double[] reflectorScalarsT, int offsetT)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dorglq(
                numRowsQ, numColsQ, numReflectors, matrixQ, offsetQ, leadingDimQ, reflectorScalarsT, offsetT,
                work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -5) info = -6;
                else if (info == -6) info = -7;
                ProcessNegativeInfo(info);
            }
        }

        internal void Dorgqr(int numRowsQ, int numColsQ, int numReflectors, double[] matrixQ, int offsetQ, int leadingDimQ,
            double[] reflectorScalarsT, int offsetT)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dorgqr(
                numRowsQ, numColsQ, numReflectors, matrixQ, offsetQ, leadingDimQ, reflectorScalarsT, offsetT,
                work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -5) info = -6;
                else if (info == -6) info = -7;
                ProcessNegativeInfo(info);
            }
        }

        internal void Dormlq(MultiplicationSide sideQ, TransposeMatrix transposeQ, int numRowsC, int numColsC, int numReflectors,
            double[] matrixQ, int offsetQ, int leadingDimQ, double[] reflectorScalarsT, int offsetT,
            double[] matrixC, int offsetC, int leadingDimC)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dormlq(
                sideQ.Translate(), transposeQ.Translate(), numRowsC, numColsC, numReflectors, matrixQ, offsetQ, leadingDimQ,
                reflectorScalarsT, offsetT, matrixC, offsetC, leadingDimC, work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -7) info = -8;
                else if (info == -8) info = -9;
                else if (info == -9) info = -11;
                else if (info == -10) info = -13;
                ProcessNegativeInfo(info);
            }
        }

        internal void Dormqr(MultiplicationSide sideQ, TransposeMatrix transposeQ, int numRowsC, int numColsC, int numReflectors,
            double[] matrixQ, int offsetQ, int leadingDimQ, double[] reflectorScalarsT, int offsetT,
            double[] matrixC, int offsetC, int leadingDimC)
        {
            int info = DefaultInfo;
            QueryWorkspaceAndExecute((work, offsetWork, lWork) => Provider.Dormqr(
                sideQ.Translate(), transposeQ.Translate(), numRowsC, numColsC, numReflectors, matrixQ, offsetQ, leadingDimQ,
                reflectorScalarsT, offsetT, matrixC, offsetC, leadingDimC, work, offsetWork, lWork, ref info));

            if (info < 0) // info can only be 0 or negative
            {
                // The indices of negative pivots must take into account the offset parameters
                if (info == -7) info = -8;
                else if (info == -8) info = -9;
                else if (info == -9) info = -11;
                else if (info == -10) info = -13;
                ProcessNegativeInfo(info);
            }
        }
    }
}
