using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.Exceptions;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    static class Preconditions
    {
        public static void CheckVectorDimensions(IVectorView vector1, IVectorView vector2)
        {
            if (vector1.Length != vector2.Length)
                throw new NonMatchingDimensionsException(string.Format(
                    "Vector1 has length of {0}, while vector2 has length of {1}", vector1.Length, vector2.Length));
        }
    }
}
