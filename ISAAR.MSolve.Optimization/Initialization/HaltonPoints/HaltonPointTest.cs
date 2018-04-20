using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Initialization.HaltonPoints
{
    public class HaltonPointTest
    {
        private static double tolerance = 1e-8;
        private static double[] testSequence1 = 
            { 1.0 / 2, 1.0 / 4, 3.0 / 4, 1.0 / 8, 5.0 / 8, 3.0 / 8, 7.0 / 8, 1.0 / 16, 9.0 / 16 };
        private static double[] testSequence2 =
            { 1.0 / 3, 2.0 / 3, 1.0 / 9, 4.0 / 9, 7.0 / 9, 2.0 / 9, 5.0 / 9, 8.0 / 9, 1.0 / 27 };

        public static void Run()
        {
            var halton = new HaltonPointGenerator(2);
            for (int i = 0; i < 9; ++i)
            {
                double[] point = halton.NextPoint();
                if (AreDifferent(point[0], testSequence1[i]) || 
                    AreDifferent(point[1], testSequence2[i]))
                {
                    Console.WriteLine("Test was unsuccessful");
                }
                //Console.WriteLine("Index = {0}: x = {1}, y = {2}", i, point[0], point[1]);
            }
            Console.WriteLine("Test was successful");
        }

        private static bool AreDifferent(double x, double y)
        {
            return (Math.Abs(1.0 - y / x) < tolerance)? false : true;
        }
    }
}
