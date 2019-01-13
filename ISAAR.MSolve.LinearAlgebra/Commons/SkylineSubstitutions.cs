using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Implementations of back and forward substitution operations with a matrix stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SkylineSubstitutions
    {
        internal static void SolveLinearSystem(int order, double[] skyValues, int[] skyDiagOffsets, 
            IVectorView rhs, IVector solution)
        {
            var x = Vector.CreateFromArray(SolveLinearSystem(order, skyValues, skyDiagOffsets, rhs));
            solution.CopyFrom(x);
        }

        internal static double[] SolveLinearSystem(int order, double[] skyValues, int[] skyDiagOffsets,
            IVectorView rhs)
        {
            // Copied from Stavroulakis code.

            // Copy the y vector
            double[] x = rhs.CopyToArray();

            // RHS vector reduction
            int n;
            for (n = 0; n < order; n++)
            {
                int KL = skyDiagOffsets[n] + 1;
                int KU = skyDiagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double C = 0;
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        C += skyValues[KK] * x[k];
                    }
                    x[n] -= C;
                }
            }

            // Back substitution
            for (n = 0; n < order; n++) x[n] /= skyValues[skyDiagOffsets[n]];

            n = order - 1;
            for (int l = 1; l < order; l++)
            {
                int KL = skyDiagOffsets[n] + 1;
                int KU = skyDiagOffsets[n + 1] - 1;
                if (KU >= KL)
                {
                    int k = n;
                    double xn = x[n];
                    for (int KK = KL; KK <= KU; KK++)
                    {
                        k--;
                        x[k] -= skyValues[KK] * xn;
                    }
                }
                n--;
            }

            return x;
        }
    }
}
