using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    // Note: scalar division should not be implemented. The user should multiply with the 1.0/scalar. That way division
    // by 0 does not need to be checked explicitly and will be caught in the client call.
    public static class ShortcutOperations
    {
        #region IVectorView extensions
        public static IVector Add(this IVectorView vector1, IVectorView vector2)
        {
            return vector1.DoPointwise(vector2, (x1, x2) => x1 + x2);
        }

        public static IVector MultiplyPointwise(this IVectorView vector1, IVectorView vector2)
        {
            return vector1.DoPointwise(vector2, (x1, x2) => x1 * x2);
        }
        
        public static IVector MultiplyScalar(this IVectorView vector1, double scalar)
        {
            return vector1.DoToAllEntries(x => x* scalar);
        }

        public static IVector Subtract(this IVectorView vector1, IVectorView vector2)
        {
            return vector1.DoPointwise(vector2, (x1, x2) => x1 - x2);
        }
        #endregion

        #region IVector extensions
        public static void AddIntoThis(this IVector vector1, IVectorView vector2)
        {
            vector1.DoPointwiseIntoThis(vector2, (x1, x2) => x1 + x2);
        }

        public static void MultiplyPointwiseIntoThis(this IVector vector1, IVectorView vector2)
        {
            vector1.DoPointwiseIntoThis(vector2, (x1, x2) => x1 * x2);
        }

        public static void MultiplyScalarIntoThis(this IVector vector, double scalar)
        {
            vector.DoToAllEntriesIntoThis(x => x * scalar);
        }

        public static void SubtractIntoThis(this IVector vector1, IVectorView vector2)
        {
            vector1.DoPointwiseIntoThis(vector2, (x1, x2) => x1 - x2);
        }
        #endregion

        #region IMatrixView extensions
        public static IMatrix Add(this IMatrixView matrix1, IMatrixView matrix2)
        {
            return matrix1.DoPointwise(matrix2, (x1, x2) => x1 + x2);
        }

        public static IMatrix MultiplyScalar(this IMatrixView matrix, double scalar)
        {
            return matrix.DoToAllEntries(x => x * scalar);
        }

        public static IMatrix Subtract(this IMatrixView matrix1, IMatrixView matrix2)
        {
            return matrix1.DoPointwise(matrix2, (x1, x2) => x1 - x2);
        }
        #endregion

        #region IMatrix extensions
        public static void AddIntoThis(this IMatrix matrix1, IMatrixView matrix2)
        {
            matrix1.DoPointwiseIntoThis(matrix2, (x1, x2) => x1 + x2);
        }

        public static void MultiplyScalarIntoThis(this IMatrix matrix, double scalar)
        {
            matrix.DoToAllEntriesIntoThis(x => x * scalar);
        }

        public static void SubtractIntoThis(this IMatrix matrix1, IMatrixView matrix2)
        {
            matrix1.DoPointwiseIntoThis(matrix2, (x1, x2) => x1 - x2);
        }
        #endregion
    }
}
