using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class GlobalMatrixChecker
    {
        private readonly string expectedMatrixPath;
        private readonly string expectedDofEnumerationPath;
        private readonly bool printIncorrectEntries;
        private readonly ValueComparer comparer;

        public GlobalMatrixChecker(string expectedMatrixPath, string expectedDofEnumerationPath, 
            double tolerance = 1e-4, bool printIncorrectEntries = true)
        {
            this.expectedMatrixPath = expectedMatrixPath;
            this.expectedDofEnumerationPath = expectedDofEnumerationPath;
            this.printIncorrectEntries = printIncorrectEntries;
            this.comparer = new ValueComparer(tolerance);
        }

        public void PrintGlobalMatrix(Model2D model, bool nodeMajorReordering = false)
        {
            Console.WriteLine("Global stiffness matrix:");
            IMatrix2D<double> globalMatrix = SingleGlobalSkylineAssembler.BuildGlobalMatrix(model);
            int[] permutation = DofReorder.OldToNewDofs(model, OutputReaders.ReadNodalDofs(expectedDofEnumerationPath));
            globalMatrix = MatrixUtilities.Reorder(globalMatrix, permutation);
            MatrixUtilities.PrintDense(globalMatrix);
        }

        public void CheckGlobalMatrix(Model2D model)
        {
            Console.WriteLine("Checking global stiffness matrix...");
            var errors = new StringBuilder("Errors at entries: ");
            bool isCorrect = true;

            // Retrieve the matrices
            Matrix2D<double> expectedMatrix = OutputReaders.ReadGlobalStiffnessMatrix(expectedMatrixPath);
            IMatrix2D<double> actualMatrix = SingleGlobalSkylineAssembler.BuildGlobalMatrix(model);
            int[] permutation = DofReorder.OldToNewDofs(model, OutputReaders.ReadNodalDofs(expectedDofEnumerationPath));
            actualMatrix = MatrixUtilities.Reorder(actualMatrix, permutation);

            // Check dimensions first
            if (actualMatrix.Rows != expectedMatrix.Rows) throw new ArgumentException("Non matching rows.");
            if (actualMatrix.Columns != expectedMatrix.Columns) throw new ArgumentException("Non matching columns.");

            // Check each entry
            for (int row = 0; row < actualMatrix.Rows; ++row)
            {
                for (int col = 0; col < actualMatrix.Columns; ++col)
                {
                    if (!comparer.AreEqual(actualMatrix[row, col], expectedMatrix[row, col]))
                    {
                        errors.Append("[").Append(row).Append(", ").Append(col).Append("] ");
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Global stiffness matrix is correct!");
            else if (printIncorrectEntries) Console.WriteLine(errors.ToString());
            else Console.WriteLine("Wrong global stiffness matrix!");

        }
    }
}
