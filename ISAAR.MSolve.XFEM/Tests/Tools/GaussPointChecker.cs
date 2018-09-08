using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class GaussPointChecker
    {
        private readonly string expectedGaussPointsPath;
        private readonly bool printIncorrectGaussPoints;
        private readonly ValueComparer valueComparer;
        private readonly GaussPoint2DComparer pointComparer;

        public GaussPointChecker(string expectedGaussPointsPath,
            double tolerance = 1e-5, bool printIncorrectGaussPoints = true)
        {
            this.expectedGaussPointsPath = expectedGaussPointsPath;
            this.printIncorrectGaussPoints = printIncorrectGaussPoints;
            this.valueComparer = new ValueComparer(tolerance);
            this.pointComparer = new GaussPoint2DComparer(tolerance);
        }

        public void PrintElementGaussPoints(IReadOnlyList<XContinuumElement2D> elements)
        {
            for (int el = 0; el < elements.Count; ++el)
            {
                //if (el == 4)
                //{
                //    Console.WriteLine("element 4: ");
                //}
                Console.WriteLine("Element " + el + ": ");

                // Retrieve the gauss points
                var actualPoints = elements[el].IntegrationStrategy.GenerateIntegrationPoints(elements[el]).ToArray();
                Array.Sort(actualPoints, pointComparer);

                for (int p = 0; p < actualPoints.Length; ++p)
                {
                    GaussPoint2D point = actualPoints[p];
                    Console.WriteLine("GP " + p + ": (" + point.Xi + " , " + point.Eta + " , " + point.Weight + ")");
                }
                Console.WriteLine();
            }
        }

        public void CheckElementGaussPoints(IReadOnlyList<XContinuumElement2D> elements)
        {
            Console.WriteLine("Checking Gauss points...");

            GaussPoint2D[][] allGaussPoints = OutputReaders.ReadAllGaussPoints(expectedGaussPointsPath, elements.Count);
            for (int el = 0; el < elements.Count; ++el)
            {
                var msg = new StringBuilder("Element " + el + ": \n");
                var errors = new StringBuilder("Different Gauss Points: \n" );
                bool isOk = true;

                // Retrieve the gauss points
                GaussPoint2D[] expectedPoints = allGaussPoints[el];
                var actualPoints = elements[el].IntegrationStrategy.GenerateIntegrationPoints(elements[el]).ToArray();

                //Check count first
                int gpCount = actualPoints.Length;
                if (expectedPoints.Length != gpCount)
                {
                    throw new Exception(msg.Append("Different number of Gauss points!").ToString());
                }

                // Sort the arrays and then check them entrywise
                Array.Sort(expectedPoints, pointComparer);
                Array.Sort(actualPoints, pointComparer);

                for (int i = 0; i < gpCount; ++i)
                {
                    GaussPoint2D p1 = expectedPoints[i];
                    GaussPoint2D p2 = actualPoints[i];
                    if ((!valueComparer.AreEqual(p1.Xi, p2.Xi)) || (!valueComparer.AreEqual(p1.Eta, p2.Eta))
                        || (!valueComparer.AreEqual(p1.Weight, p2.Weight)))
                    {
                        isOk = false;
                        if (printIncorrectGaussPoints)
                        {
                            errors.Append("GP " + i + ":");
                            errors.Append("\texpected(" + p1.Xi + " , " + p1.Eta + " , " + p1.Weight + ")\n");
                            errors.Append("        actual(" + p2.Xi + " , " + p2.Eta + " , " + p2.Weight + ")");
                            errors.Append("\n");
                        }
                        else  errors.Append(i + " ");
                    }
                }

                if (isOk) msg.Append("OK");
                else msg.Append(errors);
                Console.WriteLine(msg.ToString());
                Console.WriteLine();
            }
        }
    }
}
