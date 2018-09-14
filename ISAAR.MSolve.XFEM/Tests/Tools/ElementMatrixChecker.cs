using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class ElementMatrixChecker
    {
        private readonly string expectedMatricesPath;
        private readonly bool printIncorrectElementMatrices;
        private readonly ValueComparer comparer;

        public ElementMatrixChecker(string expectedMatricesPath, 
            double tolerance = 1e-4, bool printIncorrectElementMatrices = true)
        {
            this.expectedMatricesPath = expectedMatricesPath;
            this.printIncorrectElementMatrices = printIncorrectElementMatrices;
            this.comparer = new ValueComparer(tolerance);
        }

        public void CheckElementMatrices(IReadOnlyList<XContinuumElement2D> elements)
        {
            Console.WriteLine("Checking element stiffness matrices...");

            Matrix[] correctMatrices = 
                OutputReaders.ReadElementStiffnessMatrices(expectedMatricesPath, elements.Count);
            for (int el = 0; el < elements.Count; ++el)
            {
                var builder = new StringBuilder("Element " + el + ": \n");
                bool isOk = true;

                // Retrieve the matrices
                Matrix correctK = correctMatrices[el];
                XContinuumElement2D element = elements[el];
                Matrix kss = element.BuildStandardStiffnessMatrix();
                (Matrix kee, Matrix kes) = element.BuildEnrichedStiffnessMatricesLower();


                // Check dimensions first
                //int stdDofsCount = element.StandardFiniteElement.DofsCount;
                //int enrDofsCount = element.CountArtificialDofs();
                if (kss.NumRows + kes.NumRows != correctK.NumRows)
                    throw new ArgumentException(builder.Append("Non matching rows.").ToString());
                if (kss.NumColumns + kee.NumColumns != correctK.NumColumns)
                    throw new ArgumentException(builder.Append("Non matching columns.").ToString());

                // Check Kss entrywise
                for (int row = 0; row < kss.NumRows; ++row)
                {
                    for (int col = 0; col < kss.NumColumns; ++col)
                    {
                        if (!comparer.AreEqual(kss[row, col], correctK[row, col]))
                        {
                            isOk = false;
                            builder.Append("Error at Kss[" + row + ", " + col + "], ");
                        }
                    }
                }
                builder.Append("\n");

                // Check Kes entrywise
                for (int row = 0; row < kes.NumRows; ++row)
                {
                    for (int col = 0; col < kes.NumColumns; ++col)
                    {
                        if (!comparer.AreEqual(kes[row, col], correctK[kss.NumRows + row, col]))
                        {
                            isOk = false;
                            builder.Append("Error at Kes[" + row + ", " + col + "], ");
                        }
                    }
                }
                builder.Append("\n");

                // Check Kee entrywise
                for (int row = 0; row < kee.NumRows; ++row)
                {
                    for (int col = 0; col < kee.NumColumns; ++col)
                    {
                        if (!comparer.AreEqual(kee[row, col], correctK[kss.NumRows + row, kss.NumColumns + col]))
                        {
                            isOk = false;
                            builder.Append("Error at Kee[" + row + ", " + col + "], ");
                        }
                    }
                }
                builder.Append("\n");

                if (isOk)
                {
                    Console.WriteLine("Element " + el + " is ok.");
                }
                else if (printIncorrectElementMatrices)
                {
                    Console.WriteLine(builder.ToString());
                    Console.Write("Expected ");
                    MatrixPrinter.PrintElementMatrix(el, correctK);
                    Console.Write("Actual ");
                    MatrixPrinter.PrintElementMatrices(el, kss, kes, kee);
                }
            }
        }
    }
}
