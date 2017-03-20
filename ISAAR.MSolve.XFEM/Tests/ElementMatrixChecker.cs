using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Tests.Tools;

namespace ISAAR.MSolve.XFEM.Tests
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

            Matrix2D<double>[] correctMatrices = 
                OutputReaders.ReadElementStiffnessMatrices(expectedMatricesPath, elements.Count);
            for (int el = 0; el < elements.Count; ++el)
            {
                var builder = new StringBuilder("Element " + el + ": \n");
                bool isOk = true;

                // Retrieve the matrices
                Matrix2D<double> correctK = correctMatrices[el];
                XContinuumElement2D element = elements[el];
                SymmetricMatrix2D<double> kss, kee;
                Matrix2D<double> kes;
                kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out kes, out kee);


                // Check dimensions first
                //int stdDofsCount = element.StandardFiniteElement.DofsCount;
                //int enrDofsCount = element.CountArtificialDofs();
                if (kss.Rows + kes.Rows != correctK.Rows)
                    throw new ArgumentException(builder.Append("Non matching rows.").ToString());
                if (kss.Columns + kee.Columns != correctK.Columns)
                    throw new ArgumentException(builder.Append("Non matching columns.").ToString());

                // Check Kss entrywise
                for (int row = 0; row < kss.Rows; ++row)
                {
                    for (int col = 0; col < kss.Columns; ++col)
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
                for (int row = 0; row < kes.Rows; ++row)
                {
                    for (int col = 0; col < kes.Columns; ++col)
                    {
                        if (!comparer.AreEqual(kes[row, col], correctK[kss.Rows + row, col]))
                        {
                            isOk = false;
                            builder.Append("Error at Kes[" + row + ", " + col + "], ");
                        }
                    }
                }
                builder.Append("\n");

                // Check Kee entrywise
                for (int row = 0; row < kee.Rows; ++row)
                {
                    for (int col = 0; col < kee.Columns; ++col)
                    {
                        if (!comparer.AreEqual(kee[row, col], correctK[kss.Rows + row, kss.Columns + col]))
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
