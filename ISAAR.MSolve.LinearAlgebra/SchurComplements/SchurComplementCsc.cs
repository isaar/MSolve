using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.SchurComplements
{
    public static class SchurComplementCsc
    {
        /// <summary>
        /// Calculates the Schur complement of M/C = S = A - B^T * inv(C) * B, where M = [A B; B^T C].
        /// This method constructs inv(C) * B one column at a time and uses that column to calculate the superdiagonal
        /// entries of the corresponding column of B^T * inv(C) * B.
        /// </summary>
        public static SymmetricMatrix CalcSchurComplementSymmetric(SymmetricMatrix A, CscMatrix B, ITriangulation inverseC)
        { //TODO: Unfortunately this cannot take advantage of MKL for CSC^T * vector.
            double[] valuesB = B.RawValues;
            int[] rowIndicesB = B.RawRowIndices;
            int[] colOffsetsB = B.RawColOffsets;
            var S = SymmetricMatrix.CreateZero(A.Order);

            for (int j = 0; j < B.NumColumns; ++j)
            {
                // column j of (inv(C) * B) = inv(C) * column j of B
                Vector colB = B.GetColumn(j);
                double[] colInvCB = inverseC.SolveLinearSystem(colB).RawData;

                // column j of (B^T * inv(C) * B) = B^T * column j of (inv(C) * B)
                // However we only need the superdiagonal part of this column. 
                // Thus we only multiply the rows i of B^T (stored as columns i of B) with i <= j. 
                for (int i = 0; i <= j; ++i)
                {
                    double dot = 0.0;
                    int colStart = colOffsetsB[i]; //inclusive
                    int colEnd = colOffsetsB[i + 1]; //exclusive
                    for (int k = colStart; k < colEnd; ++k) dot += valuesB[k] * colInvCB[rowIndicesB[k]];

                    // Perform the subtraction S = A - (B^T * inv(C) * B) for the current (i, j)
                    int indexS = S.Find1DIndex(i, j);
                    S.RawData[indexS] = A.RawData[indexS] - dot;
                }
            }

            return S;
        }

        public static Matrix CalcSchurComplementFull(Matrix A, CscMatrix B, LdlSkyline inverseC)
        {
            // S = A - B^T * inv(C) * B
            Matrix invCB = Matrix.CreateZero(inverseC.Order, B.NumColumns);
            inverseC.SolveLinearSystems(B, invCB);
            return A - B.MultiplyRight(invCB, true);
        }
    }
}
