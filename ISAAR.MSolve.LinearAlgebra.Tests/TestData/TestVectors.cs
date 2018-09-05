using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// Hardcoded vectors and the results of operations between them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class TestVectors
    {
        internal static readonly double[] vector1 = {
            4.602954, 4.717625, 0.344063, 3.068357, 1.371406, 4.356552, 1.298073, 1.672528, 2.855245, 0.056439 };

        internal static readonly double[] vector2 = {
            0.690595, 0.091956, 0.581702, 0.888653, 1.793449, 1.691994, 1.115221, 0.355421, 0.647626, 1.959777 };

        internal static readonly double[] sum = {
            5.293549000000001, 4.809581000000000, 0.925765000000000, 3.957010000000000, 3.164855000000000, 6.048546000000000, 2.413294000000000, 2.027949000000000, 3.502871000000000, 2.016216000000000 };

        internal static readonly double[] difference = {
            3.912359000000000, 4.625669000000000, -0.237639000000000, 2.179704000000000, -0.422043000000000, 2.664558000000000, 0.182852000000000, 1.317107000000000, 2.207619000000000, -1.903338000000000 };


        internal static readonly double[] vector1Times2 = {
            9.205908000000001, 9.435250000000000, 0.688126000000000, 6.136714000000000, 2.742812000000000, 8.713104000000000, 2.596146000000000, 3.345056000000000, 5.710490000000000, 0.112878000000000 };

        internal static readonly double scalar2 = 3.0;

        internal static readonly double[] vector1PlusVector2Times3 = {
            6.67473900000000, 4.99349300000000, 2.08916900000000, 5.73431600000000, 6.75175300000000, 9.43253400000000, 4.64373600000000, 2.73879100000000, 4.79812300000000, 5.93577000000000 };

        internal static readonly double[] hadamardProduct = {
            3.178777017630000, 0.433813924500000, 0.200142135226000, 2.726704653121000, 2.459546719294000, 7.371259844688000, 1.447638269133000, 0.594451574288000, 1.849130898370000, 0.110607854103000 };

        internal static readonly double dotProduct = 20.3720728903530;

        internal static readonly double norm2OfVector1 = 9.29917295970765;

        internal static readonly string filePath = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
           + @"\Resources\DenseVector10.txt";
    }
}
