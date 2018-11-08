//using System;
//using System.Collections.Generic;
//using System.Text;
//using ISAAR.MSolve.LinearAlgebra.Matrices;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
//using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

////TODO: Remove this as soon as the legacy matrices are purged.
//namespace ISAAR.MSolve.Solvers.Commons
//{
//    internal static class Conversions
//    {
//        internal static Matrix MatrixOldToNew(IMatrix2D oldMatrix)
//        {
//            if (oldMatrix is Matrix2D casted) return Matrix.CreateFromArray(casted.Data);
//            else
//            {
//                int m = oldMatrix.Rows;
//                int n = oldMatrix.Columns;
//                var newMatrix = Matrix.CreateZero(m, n);
//                for (int i = 0; i < m; ++i)
//                {
//                    for (int j = 0; j < n; ++j)
//                    {
//                        newMatrix[i, j] = oldMatrix[i, j];
//                    }
//                }
//                return newMatrix;
//            }
//        }

//        internal static Matrix2D MatrixNewToOld(Matrix newMatrix) => new Matrix2D(newMatrix.CopyToArray2D());
//    }
//}
