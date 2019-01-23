using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Utilities
{
    /// <summary>
    /// Utility methods for generating random matrices and vectors.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class RandomUtilities
    {
        internal static DokSymmetric CreateRandomMatrix(int order, double nonZeroChance)
        {
            var rand = new Random();
            var dok = DokSymmetric.CreateEmpty(order);
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (rand.NextDouble() <= nonZeroChance)
                    {
                        dok[i, j] = rand.NextDouble();
                    }
                }
            }
            return dok;
        }

        internal static DokRowMajor CreateRandomSparseMatrix(int numRows, int numCols, double nonZeroChance)
        {
            var rand = new Random();
            var dok = DokRowMajor.CreateEmpty(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                for (int j = 0; j < numCols; ++j)
                {
                    if (rand.NextDouble() <= nonZeroChance)
                    {
                        dok[i, j] = rand.NextDouble();
                    }
                }
            }
            return dok;
        }

        internal static Vector CreateRandomVector(int length)
        {
            var rand = new Random();
            var vector = new double[length];
            for (int i = 0; i < length; ++i)
            {
                vector[i] = rand.NextDouble();
            }
            return Vector.CreateFromArray(vector, false);
        }
    }
}
