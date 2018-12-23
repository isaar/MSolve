using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// Hardcoded vectors and the results of operations between them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class TestVectors
    {
        internal static double[] Vector1 => new double[] {
            4.602954, 4.717625, 0.344063, 3.068357, 1.371406, 4.356552, 1.298073, 1.672528, 2.855245, 0.056439 };

        internal static double[] Vector2 => new double[] {
            0.690595, 0.091956, 0.581702, 0.888653, 1.793449, 1.691994, 1.115221, 0.355421, 0.647626, 1.959777 };

        internal static double[] Sum => new double[] {
            5.293549, 4.809581, 0.925765, 3.957010, 3.164855, 6.048546, 2.413294, 2.027949, 3.502871, 2.016216 };

        internal static double[] Difference => new double[] {
            3.912359, 4.625669, -0.237639, 2.179704, -0.422043, 2.664558, 0.182852, 1.317107, 2.207619, -1.903338 };

        internal static double[] Vector1Times2 => new double[] {
            9.205908, 9.43525, 0.688126, 6.136714, 2.742812, 8.713104, 2.596146, 3.345056, 5.71049, 0.112878 };

        internal const double Scalar2 = 3.0;

        internal static double[] Vector1PlusVector2Times3 => new double[] {
            6.674739, 4.993493, 2.089169, 5.734316, 6.751753, 9.432534, 4.643736, 2.738791, 4.798123, 5.93577 };

        internal static double[] HadamardProduct => new double[] {
            3.17877701763, 0.4338139245, 0.200142135226, 2.726704653121, 2.459546719294, 7.371259844688, 1.447638269133, 0.594451574288, 1.84913089837, 0.110607854103 };

        internal const double DotProduct = 20.3720728903530;

        internal const double Norm2OfVector1 = 9.29917295970765;

        internal static string FilePath => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\DenseVector10.txt";
    }
}
