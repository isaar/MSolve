using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;

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
            // Read the data
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

            // Calculate expected solution
            Vector expectedSolution;
            using (CholeskySuiteSparse factorization = expectedK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                expectedSolution = factorization.SolveLinearSystem(rhs);
            }

            // Build the previous K
            int order = expectedK.NumColumns;
            var previousK = DOKSymmetricColMajor.CreateIdentity(order);
            foreach (var (row, col, val) in expectedK.EnumerateNonZerosSuperDiagonal()) previousK[row, col] = val;
            foreach (int row in addedRows) previousK.SetColumnToIdentity(row);
            //foreach (int row in removedRows) previousK.SetColumnToIdentity(row); // not sure
            //TODO: it would be better to track and recreate the previous tip rows that have been removed from current K

            // Factorize previous K, update it and calculate the solution
            Vector computedSolution;
            using (CholeskySuiteSparse factorization = previousK.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                foreach (int row in addedRows)
                {
                    SparseVector newRow = expectedK.BuildColumn(row);
                    factorization.AddRow(row, newRow);
                }
                //foreach (int row in removedRows)
                //{
                //    SparseVector newRow = expectedK.BuildColumn(row);
                //    factorization.AddRow(row, newRow);
                //}
                computedSolution = factorization.SolveLinearSystem(rhs);
            }

            // Check the solution
            double error = (computedSolution - expectedSolution).Norm2() / expectedSolution.Norm2();
            Console.WriteLine($"Normalized error = {error}");

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
