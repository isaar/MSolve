using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public static class MatrixExtensions
    {
        /// <summary>
        /// result = matrix1 + matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static Matrix Add(this Matrix matrix1, Matrix matrix2)
        {
            return matrix1.Axpy(1.0, matrix2);
        }

        /// <summary>
        /// matrix1 = matrix1 + matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void AddIntoThis(this Matrix matrix1, Matrix matrix2)
        {
            matrix1.AxpyIntoThis(1.0, matrix2);
        }

        /// <summary>
        /// result = matrix1 - matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        /// <returns></returns>
        public static Matrix Subtract(this Matrix matrix1, Matrix matrix2)
        {
            return matrix1.Axpy(-1.0, matrix2);
        }

        /// <summary>
        /// matrix1 = matrix1 - matrix2
        /// </summary>
        /// <param name="matrix1"></param>
        /// <param name="matrix2"></param>
        public static void SubtractIntoThis(this Matrix matrix1, Matrix matrix2)
        {
            matrix1.AxpyIntoThis(-1.0, matrix2);
        }

        #region linear combinations
        // Telescopic methods are much, much more convenient than using lists of matrices and coefficients
        // TODO: Perhaps I should have an axpy method (implemented using BLAS lvl1) on IMatrixView and make LinearCombination 
        // use that one, while being generic. 
        // TODO: I cannot do the same with sparse matrices, since DoEntrywise is not a closed operation. Since these will be used
        // for global matrices, I need the combination to write into one of the provided matrices.
        // TODO: Also linear combinations with other matrix types may be useful, e.g. Skyline (K) with diagonal (M), but I think 
        // that for global matrices, this should be done through concrete class to use DoEntrywiseIntoThis methods. 
        // TODO: All in all I am against any LinearCombination method. It is not that difficult to call axpy twice in a dynamic 
        // analyzer.

        public static IMatrixView LinearCombination(double coefficient1, IMatrixView matrix1, double coefficient2, 
            IMatrixView matrix2, double coefficient3, IMatrixView matrix3)
        {
            IMatrixView result = matrix1.DoEntrywise(matrix2, (x1, x2) => coefficient1 * x1 + coefficient2 * x2);
            return result.DoEntrywise(matrix3, (x, x3) => x + coefficient3 * x3);
        }

        public static IMatrixView LinearCombination(double coefficient1, IMatrixView matrix1, double coefficient2, 
            IMatrixView matrix2, double coefficient3, IMatrixView matrix3, double coefficient4, IMatrixView matrix4)
        {
            IMatrixView result = LinearCombination(coefficient1, matrix1, coefficient2, matrix2, coefficient3, matrix3);
            return result.DoEntrywise(matrix4, (x, x4) => x + coefficient4 * x4);
        }

        public static Matrix LinearCombination(double coefficient1, Matrix matrix1, double coefficient2, Matrix matrix2, 
            double coefficient3, Matrix matrix3)
        {
            Matrix result = matrix1.LinearCombination(coefficient1, coefficient2, matrix2);
            result.AxpyIntoThis(coefficient3, matrix3);
            return result;
        }

        public static Matrix LinearCombination(double coefficient1, Matrix matrix1, double coefficient2, Matrix matrix2,
            double coefficient3, Matrix matrix3, double coefficient4, Matrix matrix4) 
        {
            Matrix result = LinearCombination(coefficient1, matrix1, coefficient2, matrix2, coefficient3, matrix3);
            result.AxpyIntoThis(coefficient4, matrix4);
            return result;
        }

        public static SymmetricMatrix LinearCombination(double coefficient1, SymmetricMatrix matrix1, double coefficient2, 
            SymmetricMatrix matrix2, double coefficient3, SymmetricMatrix matrix3)
        {
            SymmetricMatrix result = matrix1.DoEntrywise(matrix2, (x1, x2) => coefficient1 * x1 + coefficient2 * x2);
            result.DoEntrywiseIntoThis(matrix3, (x, x3) => x + coefficient3 * x3);
            return result;
        }
        #endregion
    }
}
