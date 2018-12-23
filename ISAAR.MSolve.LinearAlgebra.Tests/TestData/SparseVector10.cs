using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// Hardcoded sparse vectors and the results of operations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class SparseVector10
    {
        internal const int Length = 10;

        internal static double[] VectorAsDense => new double[] {
            4.602954, 0.0, 0.0, 3.068357, 1.371406, 4.356552, 0.0, 0.0, 0.0, 0.056439 };

        internal static double[] NonZeroValues => new double[] { 4.602954, 3.068357, 1.371406, 4.356552, 0.056439 };
        internal static int[] NonZeroIndices => new int[] { 0, 3, 4, 5, 9 };

        internal static double[] OtherVector => new double[] {
            0.690595, 0.091956, 0.581702, 0.888653, 1.793449, 1.691994, 1.115221, 0.355421, 0.647626, 1.959777 };

        internal static double[] OtherVectorPlusThisVectorTimes3 => new double[] {
            14.499457, 0.091956, 0.581702, 10.093724, 5.907666999999999, 14.76165, 1.115221, 0.355421, 0.647626, 2.129094};

        internal const double DotThisTimesOther = 15.846896088836;

        internal const double Norm2OfThis = 7.173944891886611;

        //internal static string FilePath => Directory.GetParent(Directory.GetCurrentDirectory()).Parent.Parent.FullName
        //   + @"\Resources\SparseVector10.txt";
    }
}
