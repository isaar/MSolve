using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra;

//TODO: Remove this as soon as the legacy matrices are purged.
namespace ISAAR.MSolve.Solvers.Commons
{
    internal static class Conversions
    {
        internal static Matrix MatrixOldToNew(Matrix2D oldMatrix) => Matrix.CreateFromArray(oldMatrix.Data);
        internal static Matrix2D MatrixNewToOld(Matrix newMatrix) => new Matrix2D(newMatrix.CopyToArray2D());
    }
}
