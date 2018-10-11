using System;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Algorithms for various linear algebra operations that do not take into account the sparsity pattern of the matrix or 
    /// vector. These serve as a default implementation when a smarter one is not readily available.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class DenseStrategies
    {
        internal static bool AreEqual(IIndexable2D matrix1, IIndexable2D matrix2, double tolerance = 1e-13)
        {
            if ((matrix1.NumRows != matrix2.NumRows) || (matrix1.NumColumns != matrix2.NumColumns)) return false;
            var comparer = new ValueComparer(tolerance);
            for (int j = 0; j < matrix1.NumColumns; ++j)
            {
                for (int i = 0; i < matrix1.NumRows; ++i)
                {
                    if (!comparer.AreEqual(matrix1[i, j], matrix2[i, j])) return false;
                }
            }
            return true;
        }

        internal static double[,] CopyToArray2D(IIndexable2D matrix)
        {
            var result = new double[matrix.NumRows, matrix.NumColumns];
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    result[i, j] = matrix[i, j];
                }
            }
            return result;
        }

        internal static Vector DoEntrywise(IVectorView vector1, IVectorView vector2,
            Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(vector1, vector2);
            var result = Vector.CreateZero(vector1.Length);
            for (int i = 0; i < vector1.Length; ++i) result[i] = binaryOperation(vector1[i], vector2[i]);
            return result;
        }

        internal static Matrix DoEntrywise(IIndexable2D matrix1, IIndexable2D matrix2, 
            Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(matrix1, matrix2);
            var result = Matrix.CreateZero(matrix1.NumRows, matrix1.NumColumns);
            for (int j = 0; j < matrix1.NumColumns; ++j)
            {
                for (int i = 0; i < matrix1.NumRows; ++i)
                {
                    result[i, j] = binaryOperation(matrix1[i, j], matrix2[i, j]);
                }
            }
            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tolerance">Can be zero</param>
        /// <returns></returns>
        internal static bool IsZero(double[] array, double tolerance)
        {
            if (tolerance == 0)
            {
                for (int i = 0; i < array.Length; ++i)
                {
                    if (array[i] != 0.0) return false;
                }
                return true;
            }
            else
            {
                for (int i = 0; i < array.Length; ++i)
                {
                    if (Math.Abs(array[i]) > tolerance) return false;
                }
                return true;
            }
        }

        internal static Vector LinearCombination(IIndexable1D vector1, double coefficient1, IIndexable1D vector2,
            double coefficient2)
        {
            Preconditions.CheckVectorDimensions(vector1, vector2);
            var result = Vector.CreateZero(vector1.Length);
            for (int i = 0; i < vector1.Length; ++i) result[i] = coefficient1 * vector1[i] + coefficient2 * vector2[i];
            return result;
        }

        internal static Matrix LinearCombination(IIndexable2D matrix1, double coefficient1, IIndexable2D matrix2, 
            double coefficient2)
        {
            Preconditions.CheckSameMatrixDimensions(matrix1, matrix2);
            var result = Matrix.CreateZero(matrix1.NumRows, matrix1.NumColumns);
            for (int j = 0; j < matrix1.NumColumns; ++j)
            {
                for (int i = 0; i < matrix1.NumRows; ++i)
                {
                    result[i, j] = coefficient1 * matrix1[i, j] + coefficient2 * matrix2[i, j];
                }
            }
            return result;
        }

        internal static Matrix Transpose(IMatrixView matrix)
        {
            var result = Matrix.CreateZero(matrix.NumColumns, matrix.NumRows);
            for (int j = 0; j < matrix.NumColumns; ++j)
            {
                for (int i = 0; i < matrix.NumRows; ++i)
                {
                    result[j, i] = matrix[i, j];
                }
            }
            return result;
        }

        internal static Matrix Multiply(IMatrixView matrix1, IMatrixView matrix2, bool transpose1, bool transpose2)
        {
            if (transpose1)
            {
                if (transpose2)
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumRows, matrix2.NumColumns);
                    var result = Matrix.CreateZero(matrix1.NumColumns, matrix2.NumRows);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumRows; ++k)
                            {
                                result[i, j] = matrix1[k, i] * matrix2[j, k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumRows, matrix2.NumRows);
                    var result = Matrix.CreateZero(matrix1.NumColumns, matrix2.NumColumns);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumRows; ++k)
                            {
                                result[i, j] = matrix1[k, i] * matrix2[k, j];
                            }
                        }
                    }
                    return result;
                }
            }
            else
            {
                if (transpose2)
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumColumns, matrix2.NumColumns);
                    var result = Matrix.CreateZero(matrix1.NumRows, matrix2.NumRows);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k < matrix1.NumColumns; ++k)
                            {
                                result[i, j] = matrix1[i, k] * matrix2[j, k];
                            }
                        }
                    }
                    return result;
                }
                else
                {
                    Preconditions.CheckMultiplicationDimensions(matrix1.NumColumns, matrix2.NumRows);
                    var result = Matrix.CreateZero(matrix1.NumRows, matrix2.NumColumns);
                    for (int i = 0; i < result.NumRows; ++i)
                    {
                        for (int j = 0; j < result.NumColumns; ++j)
                        {
                            for (int k = 0; k <matrix1.NumColumns; ++k)
                            {
                                result[i, j] = matrix1[i, k] * matrix2[k, j];
                            }
                        }
                    }
                    return result;
                }
            }
        }

        internal static Vector Multiply(IMatrixView matrix, IVectorView vector, bool transposeMatrix)
        {
            if (transposeMatrix)
            {
                Preconditions.CheckMultiplicationDimensions(matrix.NumRows, vector.Length);
                var result = Vector.CreateZero(matrix.NumColumns);
                for (int i = 0; i < result.Length; ++i)
                {
                    for (int j = 0; j < vector.Length; ++j)
                    {
                        result[i] = matrix[j, i] * vector[j];
                    }
                }
                return result;
            }
            else
            {
                Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, vector.Length);
                var result = Vector.CreateZero(matrix.NumRows);
                for (int i = 0; i < result.Length; ++i)
                {
                    for (int j = 0; j < vector.Length; ++j)
                    {
                        result[i] = matrix[j, i] * vector[j];
                    }
                }
                return result;
            }
        }


        // This will be replaced with MKL SVD. 
        internal static void SVD(Matrix matrix, double[] w, double[,] v)
        {
            //      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
            //      double precision dsqrt,dmax1,dabs,dsign

            int i, j, k, l, ii, i1, kk, k1, ll, l1, mn, its;
            double c, f, g, h, s, x, y, z, scale, anorm, myTemp;
            bool bEsc1;

            int m = matrix.NumRows;
            int n = matrix.NumColumns;
            double[,] u = matrix.CopyToArray2D();
            //double[,] u = new double[m, n];
            //double[,] a = data as double[,];
            //Array.Copy(a, u, a.GetLength(0) * a.GetLength(1));

            bool matu = false;
            bool matv = true;
            double[] rv1 = new double[n];
            //  ierr = 0;
            //     .......... householder reduction to bidiagonal form ..........
            g = 0;
            scale = 0;
            anorm = 0;
            l = 0;
            for (i = 1; i <= n; i++)
            {
                l = i + 1;
                rv1[i - 1] = scale * g;
                g = 0;
                s = 0;
                scale = 0;

                if (i <= m)
                {
                    for (k = i; k <= m; k++)
                        scale = scale + Math.Abs(u[k - 1, i - 1]);

                    if (scale != 0)
                    {
                        for (k = i; k <= m; k++)
                        {
                            u[k - 1, i - 1] = u[k - 1, i - 1] / scale;
                            s = s + Math.Pow(u[k - 1, i - 1], 2);
                        }

                        f = u[i - 1, i - 1];
                        g = -(Math.Sqrt(s));
                        if (f < 0)
                            g = g * (-1);

                        h = f * g - s;
                        u[i - 1, i - 1] = f - g;

                        if (i != n)
                        {
                            for (j = l; j <= n; j++)
                            {
                                s = 0;
                                for (k = i; k <= m; k++)
                                    s = s + u[k - 1, i - 1] * u[k - 1, j - 1];
                                f = s / h;
                                for (k = i; k <= m; k++)
                                    u[k - 1, j - 1] = u[k - 1, j - 1] + f * u[k - 1, i - 1];
                            }
                        }

                        for (k = i; k <= m; k++)
                            u[k - 1, i - 1] = scale * u[k - 1, i - 1];
                    }
                }

                w[i - 1] = scale * g;
                g = 0;
                s = 0;
                scale = 0;

                if (!((i > m) || (i == n)))
                {
                    for (k = l; k <= n; k++)
                        scale = scale + Math.Abs(u[i - 1, k - 1]);

                    if (scale != 0)
                    {
                        for (k = l; k <= n; k++)
                        {
                            u[i - 1, k - 1] = u[i - 1, k - 1] / scale;
                            s = s + Math.Pow(u[i - 1, k - 1], 2);
                        }

                        f = u[i - 1, l - 1];
                        g = -(Math.Sqrt(s));
                        if (f < 0)
                            g = g * (-1);
                        h = f * g - s;
                        u[i - 1, l - 1] = f - g;

                        for (k = l; k <= n; k++)
                            rv1[k - 1] = u[i - 1, k - 1] / h;

                        if (i != m)
                        {
                            for (j = l; j <= m; j++)
                            {
                                s = 0;
                                for (k = l; k <= n; k++)
                                    s = s + u[j - 1, k - 1] * u[i - 1, k - 1];
                                for (k = l; k <= n; k++)
                                    u[j - 1, k - 1] = u[j - 1, k - 1] + s * rv1[k - 1];
                            }
                        }

                        for (k = l; k <= n; k++)
                            u[i - 1, k - 1] = scale * u[i - 1, k - 1];
                    }
                }

                myTemp = Math.Abs(w[i - 1]) + Math.Abs(rv1[i - 1]);
                if (anorm < myTemp)
                    anorm = myTemp;
            }

            //     .......... accumulation of right-hand transformations ..........
            if (matv)
            {
                //     .......... for i=n step -1 until 1 do -- ..........
                for (ii = 1; ii <= n; ii++)
                {
                    i = n + 1 - ii;

                    if (i != n)
                    {
                        if (g != 0)
                        {
                            for (j = l; j <= n; j++)
                                //     .......... double division avoids possible underflow ..........
                                v[j - 1, i - 1] = (u[i - 1, j - 1] / u[i - 1, l - 1]) / g;

                            for (j = l; j <= n; j++)
                            {
                                s = 0;
                                for (k = l; k <= n; k++)
                                    s = s + u[i - 1, k - 1] * v[k - 1, j - 1];

                                for (k = l; k <= n; k++)
                                    v[k - 1, j - 1] = v[k - 1, j - 1] + s * v[k - 1, i - 1];
                            }
                        }

                        for (j = l; j <= n; j++)
                        {
                            v[i - 1, j - 1] = 0;
                            v[j - 1, i - 1] = 0;
                        }
                    }

                    v[i - 1, i - 1] = 1;
                    g = rv1[i - 1];
                    l = i;
                }
            }

            //     .......... accumulation of left-hand transformations ..........
            if (matu)
            {
                //     ..........for i=min(m,n) step -1 until 1 do -- ..........
                mn = n;
                if (m < n) mn = m;

                for (ii = 1; ii <= mn; ii++)
                {
                    i = mn + 1 - ii;
                    l = i + 1;
                    g = w[i - 1];

                    if (i != n)
                        for (j = l; j <= n; j++)
                            u[i - 1, j - 1] = 0;

                    if (g != 0)
                    {
                        if (i != mn)
                        {
                            for (j = l; j <= n; j++)
                            {
                                s = 0;

                                for (k = l; k <= m; k++)
                                    s = s + u[k - 1, i - 1] * u[k - 1, j - 1];
                                //     .......... double division avoids possible underflow ..........
                                f = (s / u[i - 1, i - 1]) / g;

                                for (k = i; k <= m; k++)
                                    u[k - 1, j - 1] = u[k - 1, j - 1] + f * u[k - 1, i - 1];
                            }
                        }

                        for (j = i; j <= m; j++)
                            u[j - 1, i - 1] = u[j - 1, i - 1] / g;
                    }
                    else
                    {
                        for (j = i; j <= m; j++)
                            u[j - 1, i - 1] = 0;
                    }

                    u[i - 1, i - 1] = u[i - 1, i - 1] + 1;
                }
            }

            //     .......... diagonalization of the bidiagonal form ..........
            //     .......... for k=n step -1 until 1 do -- ..........
            for (kk = 1; kk <= n; kk++)
            {
                k1 = n - kk;
                k = k1 + 1;
                its = 0;

                //     .......... test for convergence ..........
                // 520    .......... test for splitting.
                //                for l=k step -1 until 1 do -- ..........
                //    while (l != k)
                l1 = 0;
                for (; ; )
                {
                    bEsc1 = false;
                    for (ll = 1; ll <= k; ll++)
                    {
                        l1 = k - ll;
                        l = l1 + 1;

                        if ((Math.Abs(rv1[l - 1]) + anorm) == anorm)
                        {
                            bEsc1 = true;
                            break;
                        }
                        //     .......... rv1(1) is always zero, so there is no exit
                        //                through the bottom of the loop ..........

                        if ((Math.Abs(w[l1 - 1]) + anorm) == anorm)
                            break;
                    }

                    if (!(bEsc1))
                    {
                        //     .......... cancellation of rv1(l) if l greater than 1 ..........
                        c = 0;
                        s = 1.0;

                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i - 1];
                            rv1[i - 1] = c * rv1[i - 1];
                            if ((Math.Abs(f) + anorm) == anorm)
                                break;

                            g = w[i - 1];
                            h = Math.Sqrt(f * f + g * g);
                            w[i - 1] = h;
                            c = g / h;
                            s = -f / h;
                            if (matu)
                                for (j = 1; j <= m; j++)
                                {
                                    y = u[j - 1, l1 - 1];
                                    z = u[j - 1, i - 1];
                                    u[j - 1, l1 - 1] = y * c + z * s;
                                    u[j - 1, i - 1] = -y * s + z * c;
                                }
                        }
                    }

                    z = w[k - 1];
                    if (l == k)
                        break;
                    //     .......... shift from bottom 2 by 2 minor ..........
                    //      if (its .eq. 30) go to 1000
                    its = its + 1;
                    x = w[l - 1];
                    y = w[k1 - 1];
                    g = rv1[k1 - 1];
                    h = rv1[k - 1];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                    g = Math.Sqrt(f * f + 1);
                    /*      myTemp = fabs(g);
                          if (f < 0)
                            myTemp = myTemp * (-1); */

                    if (f < 0)
                        myTemp = -g;
                    else
                        myTemp = g;
                    f = ((x - z) * (x + z) + h * (y / (f + myTemp) - h)) / x;
                    //     .......... next qr transformation ..........
                    c = 1;
                    s = 1;

                    for (i1 = l; i1 <= k1; i1++)
                    {
                        i = i1 + 1;
                        g = rv1[i - 1];
                        y = w[i - 1];
                        h = s * g;
                        g = c * g;
                        z = Math.Sqrt(f * f + h * h);
                        rv1[i1 - 1] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = -x * s + g * c;
                        h = y * s;
                        y = y * c;
                        if (matv)
                        {
                            for (j = 1; j <= n; j++)
                            {
                                x = v[j - 1, i1 - 1];
                                z = v[j - 1, i - 1];
                                v[j - 1, i1 - 1] = x * c + z * s;
                                v[j - 1, i - 1] = -x * s + z * c;
                            }
                        }

                        z = Math.Sqrt(f * f + h * h);
                        w[i1 - 1] = z;
                        //     .......... rotation can be arbitrary if z is zero ..........
                        if (z != 0)
                        {
                            c = f / z;
                            s = h / z;
                        }
                        f = c * g + s * y;
                        x = -s * g + c * y;
                        if (matu)
                        {
                            for (j = 1; j <= m; j++)
                            {
                                y = u[j - 1, i1 - 1];
                                z = u[j - 1, i - 1];
                                u[j - 1, i1 - 1] = y * c + z * s;
                                u[j - 1, i - 1] = -y * s + z * c;
                            }
                        }
                    }

                    rv1[l - 1] = 0;
                    rv1[k - 1] = f;
                    w[k - 1] = x;
                    //      go to 520
                }  // end while

                //     .......... convergence ..........
                if (z < 0)
                {
                    //     .......... w(k) is made non-negative ..........
                    w[k - 1] = -z;
                    if (matv)
                    {
                        for (j = 1; j <= n; j++)
                            v[j - 1, k - 1] = -v[j - 1, k - 1];
                    }
                }
            }
        }
    }
}
