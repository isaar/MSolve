using System.Collections.Generic;

namespace ISAAR.MSolve.LinearAlgebra.Tests.TestData
{
    /// <summary>
    /// Example of adding submatrices to a large matrix (e.g. like in FEM)
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class GlobalMatrixAssembly
    {
        internal const int subOrder = 4;
        internal const int globalOrder = 8;

        internal static readonly double[,] subMatrix1 =
        {
            { 20.1,  1.1,  2.1,  3.1 },
            {  1.1, 20.2,  2.2,  3.2 },
            {  2.1,  2.2, 20.3,  3.3 },
            {  3.1,  3.2,  3.3, 20.4 }
        };

        internal static readonly int[] globalIndices1 = { 0, 1, 2, 3 };

        internal static readonly Dictionary<int, int> indicesDictionary1 = new Dictionary<int, int> {
            { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 3 } };

        internal static readonly double[,] subMatrix2 =
        {
            { 30.1,  1.1,  2.1,  3.1 },
            {  1.1, 30.2,  2.2,  3.2 },
            {  2.1,  2.2, 30.3,  3.3 },
            {  3.1,  3.2,  3.3, 30.4 }
        };

        internal static readonly int[] globalIndices2 = { 2, 3, 4, 5 };

        internal static readonly Dictionary<int, int> indicesDictionary2 = new Dictionary<int, int> {
            { 0, 2 }, { 1, 3 }, { 2, 4 }, { 3, 5 } };

        internal static readonly double[,] subMatrix3 =
        {
            { 40.1,  1.1,  2.1,  3.1 },
            {  1.1, 40.2,  2.2,  3.2 },
            {  2.1,  2.2, 40.3,  3.3 },
            {  3.1,  3.2,  3.3, 40.4 }
        };
       
        internal static readonly int[] globalIndices3 = { 4, 5, 6, 7 };

        internal static readonly Dictionary<int, int> indicesDictionary3 = new Dictionary<int, int> {
            { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };

        internal static readonly double[,] globalMatrix =
        {
            { 20.1,  1.1,  2.1,  3.1,  0.0,  0.0,  0.0,  0.0 },
            {  1.1, 20.2,  2.2,  3.2,  0.0,  0.0,  0.0,  0.0 },
            {  2.1,  2.2, 50.4,  4.4,  2.1,  3.1,  0.0,  0.0 },
            {  3.1,  3.2,  4.4, 50.6,  2.2,  3.2,  0.0,  0.0 },
            {  0.0,  0.0,  2.1,  2.2, 70.4,  4.4,  2.1,  3.1 },
            {  0.0,  0.0,  3.1,  3.2,  4.4, 70.6,  2.2,  3.2 },
            {  0.0,  0.0,  0.0,  0.0,  2.1,  2.2, 40.3,  3.3 },
            {  0.0,  0.0,  0.0,  0.0,  3.1,  3.2,  3.3, 40.4 }
        };
    }
}
