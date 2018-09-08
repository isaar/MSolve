using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class SubmatrixChecker
    {
        private readonly bool printIncorrectEntries;
        private readonly ValueComparer comparer;

        public SubmatrixChecker(double tolerance = 1e-4, bool printIncorrectEntries = true)
        {
            this.printIncorrectEntries = printIncorrectEntries;
            this.comparer = new ValueComparer(tolerance);
        }

        public void Check(IIndexable2D expected, IIndexable2D actual)
        {
            var errors = new StringBuilder("Errors at entries: ");
            bool isCorrect = true;

            if (actual.NumRows != expected.NumRows)
                throw new ArgumentException("The 2 matrices have non matching rows.");
            if (actual.NumColumns != expected.NumColumns)
                throw new ArgumentException("The 2 matrices have non matching columns.");

            for (int row = 0; row < actual.NumRows; ++row)
            {
                for (int col = 0; col < actual.NumColumns; ++col)
                {
                    if (!comparer.AreEqual(actual[row, col], expected[row, col]))
                    {
                        errors.Append("[").Append(row).Append(", ").Append(col).Append("] ");
                        isCorrect = false;
                    }
                }
            }

            if (isCorrect) Console.WriteLine("The 2 matrices are equal!");
            else if (printIncorrectEntries) Console.WriteLine(errors.Append("\n").ToString());
            else Console.WriteLine("The 2 matrices are NOT equal!");
        }

        public void Check(IReadOnlyList<double> expected, IReadOnlyList<double> actual)
        {
            var errors = new StringBuilder("Errors at entries: ");
            bool isCorrect = true;

            if (actual.Count != expected.Count)
                throw new ArgumentException("The 2 lists have non matching size.");

            for (int i = 0; i < actual.Count; ++i)
            {
                if (!comparer.AreEqual(actual[i], expected[i]))
                {
                    errors.Append("[").Append(i).Append("] ");
                    isCorrect = false;
                }
            }

            if (isCorrect) Console.WriteLine("The 2 lists are equal!");
            else if (printIncorrectEntries) Console.WriteLine(errors.Append("\n").ToString());
            else Console.WriteLine("The 2 lists are NOT equal!");
        }

        public Matrix ExtractRelevant(IIndexable2D fullMatrix, IReadOnlyList<int> relevantDofs) //TODO: Add a Sliceable interface for most matrices
        {
            Matrix result = Matrix.CreateZero(relevantDofs.Count, relevantDofs.Count);
            for (int r = 0; r < relevantDofs.Count; ++r)
            {
                for (int c = 0; c < relevantDofs.Count; ++c)
                {
                    result[r, c] = fullMatrix[relevantDofs[r], relevantDofs[c]];
                }
            }
            return result;
        }

        public Matrix ExtractRelevant(Matrix Kss, Matrix Kes, Matrix Kee, 
            IReadOnlyList<int> relevantStdDofs, IReadOnlyList<int> relevantEnrDofs)
        {
            int totalDofs = relevantStdDofs.Count + relevantEnrDofs.Count;
            int enrDofStart = relevantStdDofs.Count;
            Matrix result = Matrix.CreateZero(totalDofs, totalDofs);

            // Kss
            for (int r = 0; r < relevantStdDofs.Count; ++r)
            {
                for (int c = 0; c < relevantStdDofs.Count; ++c)
                {
                    result[r, c] = Kss[relevantStdDofs[r], relevantStdDofs[c]];
                }
            }

            // Kee
            for (int r = 0; r < relevantEnrDofs.Count; ++r)
            {
                for (int c = 0; c < relevantEnrDofs.Count; ++c)
                {
                    result[enrDofStart + r, enrDofStart + c] = Kee[relevantEnrDofs[r], relevantEnrDofs[c]];
                }
            }

            // Kes
            for (int r = 0; r < relevantEnrDofs.Count; ++r)
            {
                for (int c = 0; c < relevantStdDofs.Count; ++c)
                {
                    result[r + enrDofStart, c] = Kes[relevantEnrDofs[r], relevantStdDofs[c]];
                    result[c, r + enrDofStart] = result[r + enrDofStart, c];
                }
            }

            return result;
        }
    }
}
