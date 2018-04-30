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
using ISAAR.MSolve.XFEM.Solvers;
using ISAAR.MSolve.XFEM.Tests.Tools;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class ReanalysisDebugging
    {
        static readonly string folderPath = @"C:\Users\Serafeim\Desktop\GRACM\Reanalysis debugging\";

        public static void Main()
        {
            TestStep(1);
        }

        private static void TestStep(int step)
        {
            /// Read the data
            string matrixPath = folderPath + "reanalysis_expected_matrix_" + step +".txt" ;
            string previousMatrixPath = folderPath + "reanalysis_expected_matrix_" + (step-1) + ".txt";
            string rhsPath = folderPath + "reanalysis_expected_rhs_" + step + ".txt";
            string removedColsPath = folderPath + "reanalysis_removed_rows_" + step + ".txt";
            string addedColsPath = folderPath + "reanalysis_added_rows_" + step + ".txt";

            var matrixReader = new CoordinateTextFileReader();
            matrixReader.ReadFromFile(matrixPath);
            DOKSymmetricColMajor expectedK = matrixReader.ToSymmetricDOK();
            matrixReader = new CoordinateTextFileReader();
            matrixReader.ReadFromFile(previousMatrixPath);
            DOKSymmetricColMajor expectedPreviousK = matrixReader.ToSymmetricDOK();
            Vector rhs = (new FullVectorReader(true)).ReadFromFile(rhsPath);
            int[] removedCols = ReadDofs(removedColsPath);
            int[] addedCols = ReadDofs(addedColsPath);

            /// Artifical example using real matrices
            //CheckManagedUpdates(expectedK, removedCols, addedCols);
            //CheckSuiteSparseColAddition(expectedK, rhs, removedCols, addedCols);
            //CheckSuiteSparseColDeletion(expectedK, rhs, removedRowCols, addedCols);
            //CheckSuiteSparseColAdditionAndDeletion(expectedK, rhs, removedCols, addedCols);

            /// Actually do the reanalysis steps
            //CheckReanalysisMatrixUpdate(expectedK, expectedPreviousK, removedCols, addedCols);
            CheckReanalysisUpdate(expectedK, expectedPreviousK, rhs, removedCols, addedCols);
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

        private static void CheckReanalysisMatrixUpdate(DOKSymmetricColMajor expectedK, DOKSymmetricColMajor expectedPreviousK,
            int[] removedCols, int[] addedCols)
        {
            var checker = new SubmatrixChecker(1e-10, true);
            int order = expectedK.NumColumns;

            /// Copy the matrix. TODO: add DOK.Copy()
            var matrix = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedPreviousK.EnumerateNonZerosSuperDiagonal()) matrix[row, col] = val;

            /// Delete rows the matrix in C#
            foreach (int col in removedCols) matrix.SetColumnToIdentity(col);

            /// Add each new col all at once. It works, but not it is different from the reanalysis process.
            //foreach (int col in addedCols) matrix.SetColumn(col, expectedK.BuildColumn(col));

            /// Add new cols incrementally, namely by ignoring the rest new cols.
            var tabooRows = new HashSet<int>(addedCols);
            //tabooRows.UnionWith(removedCols); // not sure about this one. Might be needed if row deletion happens after addition
            foreach (int col in addedCols)
            {
                tabooRows.Remove(col);
                SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, tabooRows);
                matrix.SetColumn(col, newColVector);
            }

            /// Check it
            Console.WriteLine("Checking the recreated matrix with the expected one.");
            checker.Check(expectedK, matrix); 
            //TODO: Fix it ASAP. Errors in entries that should not have been affected by the update. Perhaps the level sets are
            //      updated incorrectly

            //Console.WriteLine("Checking the DOK after deleting and adding the required rows.");
            //foreach (int row in addedCols)
            //{
            //    CheckCol(currentK, row, expectedK.BuildColumn(row));
            //}
        }

        private static void CheckReanalysisUpdate(DOKSymmetricColMajor expectedK, DOKSymmetricColMajor expectedPreviousK,
            Vector rhs, int[] removedCols, int[] addedCols)
        {
            /// Calculate expected solution
            Vector expectedSolution;
            using (CholeskySuiteSparse factorization = expectedK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Factorize previous K and update it
            Vector cholmodSolution;
            using (CholeskySuiteSparse factorization = expectedPreviousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                var tabooRows = new HashSet<int>(addedCols);
                tabooRows.UnionWith(removedCols); // not sure about this one. Might be needed if row deletion happens after addition

                foreach (int col in removedCols) factorization.DeleteRow(col);

                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    SparseVector newRowVector = ReanalysisSolver.BuildNewCol(expectedK, col, tabooRows);
                    factorization.AddRow(col, newRowVector);
                }
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Check the solution
             double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking reanalysis update: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseRowColAddition(DOKSymmetricColMajor expectedK, Vector rhs, 
            int[] removedCols, int[] addedCols)
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
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);

            /// Factorize previous K, update it and check the solution
            Vector cholmodSolution;
            using (CholeskySuiteSparse factorization = previousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                var remainingRows = new HashSet<int>(addedCols);
                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    remainingRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    SparseVector newRowVector = ReanalysisSolver.BuildNewCol(expectedK, col, remainingRows);
                    factorization.AddRow(col, newRowVector);
                }
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking col addition: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseColAdditionAndDeletion(DOKSymmetricColMajor expectedK, Vector rhs, 
            int[] removedCols, int[] addedCols)
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
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);
            foreach (int col in removedCols) previousK.SetColumnToIdentity(col); // not sure

            ///TODO: it would be better to track and recreate the previous tip rows that have been removed from current K

            /// Factorize previous K, update it and check the solution
            Vector cholmodSolution;
            using (CholeskySuiteSparse factorization = previousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                foreach (int row in removedCols.OrderBy(row => row))
                {
                    factorization.DeleteRow(row);
                }

                var remainingCols = new HashSet<int>(addedCols);
                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    remainingCols.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, remainingCols);
                    factorization.AddRow(col, newColVector);
                }
                //foreach (int row in removedRows)
                //{
                //    SparseVector newRow = expectedK.BuildColumn(row);
                //    factorization.AddRow(row, newRow);
                //}
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking col deletion & addition: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseColDeletion(DOKSymmetricColMajor expectedK, Vector rhs, 
            int[] removedCols, int[] addedCols)
        {
            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) previousK[row, col] = val;
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);

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
                foreach (int col in addedCols.OrderBy(row => row))
                {
                    factorization.DeleteRow(col);
                }
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking col deletion: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckRow(IIndexable2D matrix, int colIdx, IVectorView wholeCol, double tolerance = 1e-10)
        {
            Preconditions.CheckIndexRow(matrix, colIdx);
            Preconditions.CheckSameColDimension(matrix, wholeCol);

            var comparer = new ValueComparer(tolerance);
            var wrongRows = new List<int>();
            var matrixEntries = new List<double>();
            var vectorEntries = new List<double>();
            for (int j = 0; j < wholeCol.Length; ++j)
            {
                if (!comparer.AreEqual(matrix[colIdx, j], wholeCol[j]))
                {
                    wrongRows.Add(j);
                    matrixEntries.Add(matrix[colIdx, j]);
                    vectorEntries.Add(wholeCol[j]);
                }
            }

            if (wrongRows.Count == 0) Console.WriteLine($"Col {colIdx} of the matrix is the same as the provided vector.");
            else
            {
                Console.WriteLine($"Col {colIdx} of the matrix is the same different than the provided vector, at entries");
                Console.WriteLine(" row   matrix entry      vector entry");
                for (int t = 0; t < wrongRows.Count; ++t)
                {
                    Console.WriteLine($"{wrongRows[t]}   {matrixEntries[t]}   {vectorEntries[t]}");
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
