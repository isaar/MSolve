using System;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

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

        public void PrintGlobalMatrix(Model2D model, IDofOrderer dofOrderer, bool nodeMajorReordering = false)
        {
            Console.WriteLine("Global stiffness matrix:");
            //SingleGlobalSkylineAssembler.BuildGlobalMatrix(model, out Kff, out Kfc);
            (DokSymmetric Kff, DokRowMajor Kfc) = (new GlobalDOKAssembler()).BuildGlobalMatrix(model, dofOrderer);
            int[] permutation = DofReorder.OldToNewDofs(model, OutputReaders.ReadNodalDofs(expectedDofEnumerationPath), dofOrderer);
            (new FullMatrixWriter()).WriteToConsole(Kff.Reorder(permutation, true));
        }

        public void CheckGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            Console.WriteLine("Checking global stiffness matrix...");
            var errors = new StringBuilder("Errors at entries: ");
            bool isCorrect = true;

            // Retrieve the matrices
            IMatrixView expectedMatrix = OutputReaders.ReadGlobalStiffnessMatrix(expectedMatrixPath);
            //SingleGlobalSkylineAssembler.BuildGlobalMatrix(model, out Kff, out Kfc);
            (DokSymmetric Kff, DokRowMajor Kfc) = (new GlobalDOKAssembler()).BuildGlobalMatrix(model, dofOrderer);
            int[] permutation = DofReorder.OldToNewDofs(model, OutputReaders.ReadNodalDofs(expectedDofEnumerationPath), dofOrderer);
            IMatrixView actualMatrix = Kff.Reorder(permutation, true);

            // Check dimensions first
            if (actualMatrix.NumRows != expectedMatrix.NumRows)
                throw new ArgumentException("The 2 global matrices have non matching rows.");
            if (actualMatrix.NumColumns != expectedMatrix.NumColumns)
                throw new ArgumentException("The 2 global matrices have non matching columns.");

            // Check each entry
            for (int row = 0; row < actualMatrix.NumRows; ++row)
            {
                for (int col = 0; col < actualMatrix.NumColumns; ++col)
                {
                    if (!comparer.AreEqual(actualMatrix[row, col], expectedMatrix[row, col]))
                    {
                        errors.Append("[").Append(row).Append(", ").Append(col).Append("] ");
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Global stiffness matrix is correct!\n");
            else if (printIncorrectEntries) Console.WriteLine(errors.Append("\n").ToString());
            else Console.WriteLine("Wrong global stiffness matrix!\n");

        }
    }
}
