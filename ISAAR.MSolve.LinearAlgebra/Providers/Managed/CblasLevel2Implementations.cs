namespace ISAAR.MSolve.LinearAlgebra.Providers.Managed
{
    /// <summary>
    /// Custom and unoptimized managed implementations of BLAS like operations, for which I have not found 3rd party 
    /// implementations yet. If performance is of concern, then the user should choose the native optimized providers (e.g. MKL).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class CblasLevel2Implementations
    {
        internal enum Diagonal { Regular, Unit, Zero};

        /// <summary>
        /// x = inv(U) * x
        /// </summary>
        internal static void BackSubstitutionPackedColMajor(bool unit, int n, double[] a, int offsetA, 
            double[] x, int offsetX, int incX)
        {
            for (int i = n-1; i >= 0; --i)
            {
                double dot = 0;
                for (int j = i+1; j < n; ++j)
                {
                    // index = i + ((j + 1) * j) / 2. They are not contigous during the dot products.
                    dot += a[offsetA + i + ((j + 1) * j) / 2] * x[offsetX + j * incX];
                }

                // b[i] will not be used any more and can be replaced with x[i]
                if (unit) x[offsetX + i * incX] -= dot;
                else
                {
                    int idxX = offsetX + i * incX;
                    x[idxX] = (x[idxX] - dot) / a[offsetA + i + ((i + 1) * i) / 2];
                }
            }
        }

        /// <summary>
        /// x = inv(L) * x
        /// </summary>
        internal static void ForwardSubstitutionPackedRowMajor(bool unit, int n, double[] a, int offsetA,
            double[] x, int offsetX, int incX)
        {
            for (int i = 0; i < n; ++i)
            {
                double dot = 0;
                int rowOffset = offsetA + ((i + 1) * i) / 2;
                for (int j = 0; j < i; ++j)
                {
                    // index = j + ((i + 1) * i) / 2. They are contigous during the dot products.
                    dot += a[rowOffset + j] * x[offsetX + j * incX];
                }

                // b[i] will not be used any more and can be replaced with x[i]
                if (unit) x[offsetX + i * incX] -= dot;
                else
                {
                    int idxX = offsetX + i * incX;
                    x[idxX] = (x[idxX] - dot) / a[rowOffset + i];
                }
            }
        }

        /// <summary>
        /// y = alpha * A * x + beta * y
        /// </summary>
        internal static void LowerTimesVectorPackedRowMajor(Diagonal diag, int n, double alpha, double[] a, int offsetA,
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
        {
            for (int i = 0; i < n; ++i)
            {
                // Row * vector without the diagonal
                double dot = 0.0;
                int rowOffset = offsetA + ((i + 1) * i) / 2;
                for (int j = 0; j < i; ++j)
                {
                    // index = j + ((i + 1) * i) / 2. They are contigous during the dot products.
                    dot += a[rowOffset + j] * x[offsetX + j * incX];
                }

                // Take into account the diagonal
                if (diag == Diagonal.Regular) dot += a[rowOffset + i] * x[offsetX + i * incX];
                else if (diag == Diagonal.Unit) dot += x[offsetX + i * incX];
                // else the diagonal is ignored altogether (we operate only on the superdiagonal part)

                // Take into account the rest and store it
                int idxY = offsetY + i * incY;
                y[idxY] = alpha * dot + beta * y[idxY];
            }
        }

        /// <summary>
        /// y = alpha * A * x + beta * y
        /// </summary>
        internal static void UpperTimesVectorPackedColMajor(Diagonal diag, int n, double alpha, double[] a, int offsetA, 
            double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY)
        {
            for (int i = 0; i < n; ++i)
            {
                // Row * vector without the diagonal
                double dot = 0.0;
                for (int j = i + 1; j < n; ++j)
                {
                    // index = i + ((j + 1) * j) / 2. They are not contigous during the dot products.
                    dot += a[i + ((j + 1) * j) / 2] * x[offsetX + j * incX];
                }

                // Take into account the diagonal
                if (diag == Diagonal.Regular) dot += a[offsetA + i + ((i + 1) * i) / 2] * x[offsetX + i * incX];
                else if (diag == Diagonal.Unit) dot += x[offsetX + i * incX];
                // else the diagonal is ignored altogether (we operate only on the superdiagonal part)

                // Take into account the rest and store it
                int idxY = offsetY + i * incY;
                y[idxY] = alpha * dot + beta * y[idxY];
            }
        }
    }
}
