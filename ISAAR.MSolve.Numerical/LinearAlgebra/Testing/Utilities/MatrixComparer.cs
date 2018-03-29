using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities
{
    
    public class MatrixComparer
    {
        private readonly ValueComparer valueComparer;

        public MatrixComparer(double tolerance = 1e-6)
        {
            this.valueComparer = new ValueComparer(tolerance);
        }

        /// <summary>
        /// Prints to console all entries that are different.
        /// </summary>
        public static void CheckSameMatrices(IIndexable2D computed, IMatrix2D expected)
        {
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.NumRows != expected.Rows) || (computed.NumColumns != expected.Columns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.NumRows; ++i)
            {
                for (int j = 0; j < computed.NumColumns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                    }
                }
            }
        }

        /// <summary>
        /// Prints to console all entries that are different.
        /// </summary>
        public static void CheckSameMatrices(IMatrix2D computed, IMatrix2D expected)
        {
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.Rows != expected.Rows) || (computed.Columns != expected.Columns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.Rows; ++i)
            {
                for (int j = 0; j < computed.Columns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                    }
                }
            }
        }
    }
}
