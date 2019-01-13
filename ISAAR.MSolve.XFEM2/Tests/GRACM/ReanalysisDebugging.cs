using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Providers;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Tests.Tools;

namespace ISAAR.MSolve.XFEM.Tests.GRACM
{
    class ReanalysisDebugging
    {
        static readonly string folderPath = @"C:\Users\Serafeim\Desktop\GRACM\Reanalysis_debugging\";

        public static void Run()
        {
            TestStep(5);
        }

        private static void TestStep(int step)
        {
            /// Read the data
            string matrixPath = folderPath + "reanalysis_expected_matrix_" + step + ".txt";
            string previousMatrixPath = folderPath + "reanalysis_expected_matrix_" + (step - 1) + ".txt";
            string rhsPath = folderPath + "reanalysis_expected_rhs_" + step + ".txt";
            string removedColsPath = folderPath + "reanalysis_removed_rows_" + step + ".txt";
            string addedColsPath = folderPath + "reanalysis_added_rows_" + step + ".txt";

            var matrixReader = new CoordinateTextFileReader();
            DokSymmetric expectedK = matrixReader.ReadFileAsDokSymmetricColMajor(matrixPath);
            DokSymmetric expectedPreviousK = matrixReader.ReadFileAsDokSymmetricColMajor(previousMatrixPath);
            var rhs = Vector.CreateFromArray((new Array1DReader(true)).ReadFile(rhsPath));
            int[] removedCols = ReadDofs(removedColsPath);
            int[] addedCols = ReadDofs(addedColsPath);

            /// Artifical example using real matrices
            CheckManagedUpdates(expectedK, expectedPreviousK, rhs, removedCols, addedCols);
            //CheckManagedDelete(expectedK, expectedPreviousK, rhs, removedCols, addedCols);
            //CheckSuiteSparseColAddition(expectedK, rhs, removedCols, addedCols);
            //CheckSuiteSparseColDeletion(expectedK, rhs, removedCols, addedCols);
            //CheckSuiteSparseColAdditionAndDeletion(expectedK, rhs, removedCols, addedCols);

            /// Actually do the reanalysis steps
            //CheckReanalysisMatrixUpdate(expectedK, expectedPreviousK, removedCols, addedCols);
            //CheckReanalysisUpdate(expectedK, expectedPreviousK, rhs, removedCols, addedCols);
        }

        /// <summary>
        /// Delete all updated columns from both the new and previous matrix, in order to check the other entries.
        /// </summary>
        private static void CheckManagedDelete(DokSymmetric expectedNewK, DokSymmetric expectedPreviousK,
            Vector rhs, int[] colsToRemove, int[] colsToAdd)
        {
            var checker = new SubmatrixChecker(1e-10, true);

            /// Delete all updated dofs from both matrices
            var previousK = DokSymmetric.CreateFromSparseMatrix(expectedPreviousK);
            var newK = DokSymmetric.CreateFromSparseMatrix(expectedNewK);
            var updatedCols = new HashSet<int>(colsToAdd);
            updatedCols.UnionWith(colsToRemove);
            foreach (int col in updatedCols)
            {
                previousK.SetColumnToIdentity(col);
                newK.SetColumnToIdentity(col);
            }

            /// Check the solution vectors
            Vector previousSolution, newSolution;
            using (var factor = CholeskySuiteSparse.Factorize(
                previousK.BuildSymmetricCscMatrix(true), true))
            {
                previousSolution = factor.SolveLinearSystem(rhs);
            }
            using (var factor = CholeskySuiteSparse.Factorize(
                newK.BuildSymmetricCscMatrix(true), true))
            {
                newSolution = factor.SolveLinearSystem(rhs);
            }
            double error = (newSolution - previousSolution).Norm2() / previousSolution.Norm2();
            //Console.WriteLine($"Norm = {newSolution.Norm2()}");
            Console.WriteLine(
                $"After deleting all updated columns from the stiffness matrices of both steps, normalized error = {error}");

            /// Check the columns after deletion
            Console.WriteLine("Checking the columns of the 2 matrices after deleting them");
            foreach (var j in updatedCols)
            {
                CheckColumn(newK, previousK, j);
            }

            ///Check the columns before deletion
            //Console.WriteLine("Checking the columns of the 2 matrices before deleting them");
            //foreach (var j in updatedCols)
            //{
            //    CheckColumn(expectedNewK, expectedPreviousK, j);
            //}

            ///Check all diagonal entries before deletion
            Console.WriteLine("Checking the diagonals of the 2 matrices after deleting the columns");
            for (int i = 0; i < newK.NumRows; ++i)
            {
                if (newK[i, i] != previousK[i, i])
                {
                    Console.WriteLine($"previousK[{i}, {i}] = {previousK[i, i]} - newK[{i}, {i}] = {newK[i, i]}");
                }
            }
        }

        private static void CheckManagedUpdates(DokSymmetric expectedK, DokSymmetric previousK, Vector rhs,
            int[] colsToRemove, int[] colsToAdd)
        {
            var checker = new SubmatrixChecker(1e-10, true);

            /// Rebuild the new K from the previous one, by using the managed row add/delete methods:
            var newK = DokSymmetric.CreateFromSparseMatrix(previousK);
            var tabooRows = new HashSet<int>(colsToAdd);
            tabooRows.UnionWith(colsToRemove);
            foreach (int col in colsToRemove)
            {
                newK.SetColumnToIdentity(col);
            }
            foreach (int col in colsToAdd)
            {
                tabooRows.Remove(col);
                newK.SetColumn(col, expectedK.GetColumnWithoutRows(col, tabooRows));
            }

            /// Check the solution vectors
            Vector expectedSolution, computedSolution;
            using (var factor = CholeskySuiteSparse.Factorize(
                expectedK.BuildSymmetricCscMatrix(true), false))
            {
                expectedSolution = factor.SolveLinearSystem(rhs);
            }
            using (var factor = CholeskySuiteSparse.Factorize(
                newK.BuildSymmetricCscMatrix(true), false))
            {
                computedSolution = factor.SolveLinearSystem(rhs);
            }
            double error = (computedSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Normalized error = {error}");


            /// Check the updated columns after update
            var updatedCols = colsToAdd.Union(colsToRemove);
            Console.WriteLine("Checking the updated columns of the 2 matrices after updating them");
            foreach (var j in updatedCols)
            {
                CheckColumn(expectedK, newK, j, false);
            }

            ///Check all columns after update
            Console.WriteLine("Checking all columns of the 2 matrices after updating them");
            for(int j = 0; j < newK.NumColumns; ++j)
            {
                CheckColumn(expectedK, newK, j, false);
            }
        }

        private static void CheckManagedUpdatesOLD(DokSymmetric expectedK, int[] removedRows, int[] addedRows)
        {
            var checker = new SubmatrixChecker(1e-10, true);

            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosUpper()) previousK[row, col] = val;
            foreach (int row in addedRows) previousK.SetColumnToIdentity(row);
            //foreach (int row in removedRows) previousK.SetColumnToIdentity(row); // not sure

            /// Check that row deletion in C# works
            Console.WriteLine("Checking the DOK after deleting the required rows.");
            foreach (int row in addedRows)
            {
                var wholeRow = SparseVector.CreateFromDictionary(order, new Dictionary<int, double>() { { row, 1.0 } });
                CheckColumn(previousK, row, wholeRow);
            }
            Console.WriteLine();

            /// Recreate the new K and check it
            var currentK = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosUpper()) currentK[row, col] = val;
            //Console.WriteLine("Checking the copied matrix with the expected one.");
            //checker.Check(expectedK, currentK);
            foreach (int row in addedRows) currentK.SetColumnToIdentity(row);
            foreach (int row in addedRows) currentK.SetColumn(row, expectedK.GetColumn(row));
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

        private static void CheckReanalysisMatrixUpdate(DokSymmetric expectedK, DokSymmetric expectedPreviousK,
            int[] removedCols, int[] addedCols)
        {
            var checker = new SubmatrixChecker(1e-10, true);
            int order = expectedK.NumColumns;

            /// Copy the matrix. TODO: add DOK.Copy()
            var matrix = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedPreviousK.EnumerateNonZerosUpper()) matrix[row, col] = val;

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
                //SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, tabooRows);
                SparseVector newColVector = expectedK.GetColumnWithoutRows(col, tabooRows);
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

        private static void CheckReanalysisUpdate(DokSymmetric expectedK, DokSymmetric expectedPreviousK,
            Vector rhs, int[] removedCols, int[] addedCols)
        {
            /// Calculate expected solution
            Vector expectedSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                expectedK.BuildSymmetricCscMatrix(true), false))
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Factorize previous K and update it
            Vector cholmodSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                expectedPreviousK.BuildSymmetricCscMatrix(true), false))
            {
                var tabooRows = new HashSet<int>(addedCols);
                tabooRows.UnionWith(removedCols); // not sure about this one. Might be needed if row deletion happens after addition

                foreach (int col in removedCols) factorization.DeleteRow(col);

                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    //SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, tabooRows);
                    SparseVector newColVector = expectedK.GetColumnWithoutRows(col, tabooRows);
                    factorization.AddRow(col, newColVector);
                }
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Check the solution
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking reanalysis update: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseColAddition(DokSymmetric expectedK, Vector rhs,
            int[] removedCols, int[] addedCols)
        {
            /// Calculate expected solution
            Vector expectedSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                expectedK.BuildSymmetricCscMatrix(true), true))
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosUpper()) previousK[row, col] = val;
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);

            /// Factorize previous K, update it and check the solution
            Vector cholmodSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                previousK.BuildSymmetricCscMatrix(true), false))
            {
                var tabooRows = new HashSet<int>(addedCols);
                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    //SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, tabooRows);
                    SparseVector newColVector = expectedK.GetColumnWithoutRows(col, tabooRows);
                    factorization.AddRow(col, newColVector);
                }
                cholmodSolution = factorization.SolveLinearSystem(rhs);
            }
            double error = (cholmodSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Checking col addition: Normalized error = {error}");
            Console.WriteLine();
        }

        private static void CheckSuiteSparseColAdditionAndDeletion(DokSymmetric expectedK, Vector rhs,
            int[] removedCols, int[] addedCols)
        {
            /// Calculate expected solution
            Vector expectedSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                expectedK.BuildSymmetricCscMatrix(true), false))
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosUpper()) previousK[row, col] = val;
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);
            foreach (int col in removedCols) previousK.SetColumnToIdentity(col); // not sure

            ///TODO: it would be better to track and recreate the previous tip rows that have been removed from current K

            /// Factorize previous K, update it and check the solution
            Vector cholmodSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                previousK.BuildSymmetricCscMatrix(true), false))
            {
                foreach (int row in removedCols.OrderBy(row => row))
                {
                    factorization.DeleteRow(row);
                }

                var tabooRows = new HashSet<int>(addedCols);
                foreach (int col in addedCols)
                //foreach (int row in addedRows.OrderBy(row => row)) //Not needed, but may be faster
                {
                    tabooRows.Remove(col); //It must be called before passing the remaining cols as taboo rows.
                    //SparseVector newColVector = ReanalysisSolver.BuildNewCol(expectedK, col, remainingCols);
                    SparseVector newColVector = expectedK.GetColumnWithoutRows(col, tabooRows);
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

        private static void CheckSuiteSparseColDeletion(DokSymmetric expectedK, Vector rhs,
            int[] removedCols, int[] addedCols)
        {
            /// Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DokSymmetric.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosUpper()) previousK[row, col] = val;
            foreach (int col in addedCols) previousK.SetColumnToIdentity(col);

            /// Calculate expected solution from previous K
            Vector expectedSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                previousK.BuildSymmetricCscMatrix(true), false))
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            /// Factorize current K, update it and check the solution
            Vector cholmodSolution;
            using (var factorization = CholeskySuiteSparse.Factorize(
                expectedK.BuildSymmetricCscMatrix(true), false))
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

        private static void CheckColumn(IIndexable2D matrix1, IIndexable2D matrix2, int colIdx, bool printIfCorrect = true, 
            double tolerance = 1e-10)
        {
            Preconditions.CheckIndexCol(matrix1, colIdx);
            Preconditions.CheckIndexCol(matrix2, colIdx);
            Preconditions.CheckSameMatrixDimensions(matrix1, matrix2);

            var comparer = new ValueComparer(tolerance);
            var wrongRows = new List<int>();
            var matrix1Entries = new List<double>();
            var matrix2Entries = new List<double>();
            for (int i = 0; i < matrix1.NumRows; ++i)
            {
                if (!comparer.AreEqual(matrix1[i, colIdx], matrix2[i, colIdx]))
                {
                    wrongRows.Add(i);
                    matrix1Entries.Add(matrix1[i, colIdx]);
                    matrix2Entries.Add(matrix2[i, colIdx]);
                }
            }

            if (wrongRows.Count == 0)
            {
                if (printIfCorrect) Console.WriteLine($"Column {colIdx} is the same in both matrices.");
            }
            else
            {
                Console.WriteLine($"Columns {colIdx} of the two matrices are different, at entries");
                Console.WriteLine(" row   matrix1      matrix2");
                for (int t = 0; t < wrongRows.Count; ++t)
                {
                    Console.WriteLine($"{wrongRows[t]}   {matrix1Entries[t]}   {matrix2Entries[t]}");
                }
            }
        }

        private static void CheckColumn(IIndexable2D matrix, int colIdx, IVectorView wholeCol, double tolerance = 1e-10)
        {
            Preconditions.CheckIndexCol(matrix, colIdx);
            Preconditions.CheckSameRowDimension(matrix, wholeCol);

            var comparer = new ValueComparer(tolerance);
            var wrongRows = new List<int>();
            var matrixEntries = new List<double>();
            var vectorEntries = new List<double>();
            for (int i = 0; i < wholeCol.Length; ++i)
            {
                if (!comparer.AreEqual(matrix[i, colIdx], wholeCol[i]))
                {
                    wrongRows.Add(i);
                    matrixEntries.Add(matrix[i, colIdx]);
                    vectorEntries.Add(wholeCol[i]);
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
