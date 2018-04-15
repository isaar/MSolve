using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestLibs
{
    class TestSuiteSparseUpdate
    {
        public static void CheckRowAddition()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            VectorMKL rhs = VectorMKL.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix as diagonal
            var matrixExpected = Matrix.CreateZero(original.NumRows, original.NumColumns);
            var dok = new DOKSymmetricColMajor(SparsePositiveDefinite.order);
            for (int j = 0; j < matrixExpected.NumColumns; ++j)
            {
                matrixExpected[j, j] = 1.0;
                dok[j, j] = 1.0;
            }
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky();

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                matrixExpected.SetSubmatrix(0, 0, original.Slice(0, i + 1, 0, i + 1));
                Console.WriteLine($"\nOnly dofs [0, {i}]");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.AddRow(i, matrixExpected.SliceRow(i));

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                VectorMKL solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                VectorMKL solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.CheckVectorEquality(solutionExpected, solutionComputed);
            }
        }

        public static void CheckRowDeletion()
        {
            FullMatrixWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 1, MaxIntegerDigits = 2 };
            FullVectorWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Array1DWriter.NumericFormat = new FixedPointFormat { NumDecimalDigits = 4, MaxIntegerDigits = 3 };
            Comparer comparer = new Comparer(Comparer.PrintMode.Always);

            Matrix original = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            VectorMKL rhs = VectorMKL.CreateFromArray(SparsePositiveDefinite.rhs);
            Console.WriteLine("Full matrix: ");
            (new FullMatrixWriter(original)).WriteToConsole();

            // Start the matrix from the original
            var matrixExpected = Matrix.CreateFromArray(SparsePositiveDefinite.matrix);
            var dok = new DOKSymmetricColMajor(SparsePositiveDefinite.order);
            for (int j = 0; j < matrixExpected.NumColumns; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (matrixExpected[i, j] != 0) dok[i, j] = matrixExpected[i, j];
                }
            }
            CholeskySuiteSparse factor = dok.BuildSymmetricCSCMatrix(true).FactorCholesky();

            for (int i = 0; i < matrixExpected.NumRows; ++i)
            {
                // Update matrix
                VectorMKL identityRow = VectorMKL.CreateZero(matrixExpected.NumColumns);
                identityRow[i] = 1.0;
                matrixExpected.SetRow(i, identityRow);
                matrixExpected.SetColumn(i, identityRow);
                Console.WriteLine($"\nOnly dofs [{i+1}, 10)");
                (new FullMatrixWriter(matrixExpected)).WriteToConsole();
                factor.DeleteRow(i);

                // Solve new linear system
                Console.WriteLine("\nCheck linear system solution");
                VectorMKL solutionExpected = matrixExpected.FactorCholesky().SolveLinearSystem(rhs);
                VectorMKL solutionComputed = factor.SolveLinearSystem(rhs);
                comparer.CheckVectorEquality(solutionExpected, solutionComputed);
            }
        }
    }
}
