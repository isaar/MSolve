using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
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
    }
}
