using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Tests.Tools;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class ReanalysisDebugging
    {
        static readonly string folderPath = @"C:\Users\Serafeim\Desktop\GRACM\Reanalysis debugging\";

        public static void Main()
        {
            TestStep();
        }

        private static void TestStep()
        {
            /// Read the data
            string matrixPath = folderPath + "reanalysis_expected_matrix.txt";
            string rhsPath = folderPath + "reanalysis_expected_rhs.txt";
            string removedRowsPath = folderPath + "reanalysis_removed_rows.txt";
            string addedRowsPath = folderPath + "reanalysis_added_rows.txt";
            var matrixReader = new CoordinateTextFileReader();
            matrixReader.ReadFromFile(matrixPath);
            DOKSymmetricColMajor expectedK = matrixReader.ToSymmetricDOK();
            Vector rhs = (new FullVectorReader(true)).ReadFromFile(rhsPath);
            int[] removedRows = ReadDofs(removedRowsPath);
            int[] addedRows = ReadDofs(addedRowsPath);


            //CheckManagedUpdates(expectedK, removedRows, addedRows);
            //CheckSuiteSparseRowAddition(expectedK, rhs, removedRows, addedRows);
            CheckSuiteSparseRowDeletion(expectedK, rhs, removedRows, addedRows);
        }

        private static void CheckManagedUpdates(DOKSymmetricColMajor expectedK, int[] removedRows, int[] addedRows)
        {
            var checker = new SubmatrixChecker(1e-10, true);

            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) previousK[row, col] = val;
            foreach (int row in addedRows) previousK.SetColumnToIdentity(row);
            //foreach (int row in removedRows) previousK.SetColumnToIdentity(row); // not sure

            /// Check that row deletion in C# works
            Console.WriteLine("Checking the DOK after deleting the required rows.");
            foreach (int row in addedRows)
            {
                var wholeRow = SparseVector.CreateFromDictionary(order, new Dictionary<int, double>() { { row, 1.0 } });
                CheckRow(previousK, row, wholeRow);
            }
            Console.WriteLine();

            /// Recreate the new K and check it
            var currentK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) currentK[row, col] = val;
            //Console.WriteLine("Checking the copied matrix with the expected one.");
            //checker.Check(expectedK, currentK);
            foreach (int row in addedRows) currentK.SetColumnToIdentity(row);
            foreach (int row in addedRows) currentK.SetColumn(row, expectedK.BuildColumn(row));
            Console.WriteLine("Checking the recreated matrix with the expected one.");
            checker.Check(expectedK, currentK);
            //if (currentK.Equals(expectedK, 1e-5)) Console.WriteLine("Getting/setting rows in C# seems to work fine. Recreated matrix is the same.");
            //else Console.WriteLine("Error in getting/setting rows in C#. Recreated matrix is not the same.");
            //Console.WriteLine("Checking the DOK after deleting and adding the required rows.");
            //foreach (int row in addedRows)
            //{
            //    CheckRow(currentK, row, expectedK.BuildColumn(row));
            //}
        }

        private static void CheckSuiteSparseRowAddition(DOKSymmetricColMajor expectedK, Vector rhs, int[] removedRows, int[] addedRows)
        {
            /// Calculate expected solution
            Vector expectedSolution;
            using (CholeskySuiteSparse factorization = expectedK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) previousK[row, col] = val;
            foreach (int row in addedRows) previousK.SetColumnToIdentity(row);
            //foreach (int row in removedRows) previousK.SetColumnToIdentity(row); // not sure

            ///TODO: it would be better to track and recreate the previous tip rows that have been removed from current K

            /// Factorize previous K, update it and check the solution
            Vector cholmodSolution;
            using (CholeskySuiteSparse factorization = previousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                foreach (int row in addedRows.OrderBy(row => row))
                {
                    SparseVector newRow = expectedK.BuildColumn(row);
                    factorization.AddRow(row, newRow);
                }
                //foreach (int row in removedRows)
                //{
                //    SparseVector newRow = expectedK.BuildColumn(row);
                //    factorization.AddRow(row, newRow);
                //}
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseRowDeletion(DOKSymmetricColMajor expectedK, Vector rhs, int[] removedRows, int[] addedRows)
        {
            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) previousK[row, col] = val;
            foreach (int row in addedRows) previousK.SetColumnToIdentity(row);
            //foreach (int row in removedRows) previousK.SetColumnToIdentity(row); // not sure
            ///TODO: it would be better to track and recreate the previous tip rows that have been removed from current K

            /// Calculate expected solution from previous K
            Vector expectedSolution;
            using (CholeskySuiteSparse factorization = previousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Factorize current K, update it and check the solution
            Vector cholmodSolution;
            using (CholeskySuiteSparse factorization = expectedK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                //foreach (int row in addedRows)
                foreach (int row in addedRows.OrderBy(row => row))
                {
                    factorization.DeleteRow(row);
                }
                //foreach (int row in removedRows)
                //{
                //    SparseVector newRow = expectedK.BuildColumn(row);
                //    factorization.AddRow(row, newRow);
                //}
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckRow(IIndexable2D matrix, int rowIdx, IVectorView wholeRow, double tolerance = 1e-10)
        {
            Preconditions.CheckIndexRow(matrix, rowIdx);
            Preconditions.CheckSameColDimension(matrix, wholeRow);

            var comparer = new ValueComparer(tolerance);
            var wrongCols = new List<int>();
            var matrixEntries = new List<double>();
            var vectorEntries = new List<double>();
            for (int j = 0; j < wholeRow.Length; ++j)
            {
                if (!comparer.AreEqual(matrix[rowIdx, j], wholeRow[j]))
                {
                    wrongCols.Add(j);
                    matrixEntries.Add(matrix[rowIdx, j]);
                    vectorEntries.Add(wholeRow[j]);
                }
            }

            if (wrongCols.Count == 0) Console.WriteLine($"Row {rowIdx} of the matrix is the same as the provided vector.");
            else
            {
                Console.WriteLine($"Row {rowIdx} of the matrix is the same different than the provided vector, at entries");
                Console.WriteLine(" col   matrix entry      vector entry");
                for (int t = 0; t < wrongCols.Count; ++t)
                {
                    Console.WriteLine($"{wrongCols[t]}   {matrixEntries[t]}   {vectorEntries[t]}");
                }
            }
        }

        private static int[] ReadDofs(string path)
        {
            using (var reader = new StreamReader(path))
            {
                // Read the list length
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                int length = Int32.Parse(firstLine);
                var dofs = new int[length];

                // Read the list entries. TODO: what if each entry is on a different line?
                string line = reader.ReadLine();
                string[] subStrings = line.Split(new char[] { }, StringSplitOptions.RemoveEmptyEntries);
                if (subStrings.Length != length) throw new IOException("Mismatch in provided entries and their declared count.");
                for (int i = 0; i < length; ++i)
                {
                    dofs[i] = Int32.Parse(subStrings[i]);
                }

                return dofs;
            }
        }

    }
}
