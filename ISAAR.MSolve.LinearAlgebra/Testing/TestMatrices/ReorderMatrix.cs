using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class ReorderMatrix
    {
        public const int order = 5;

        public static readonly double[,] matrix = {
            { 11.0,  0.0,  0.0,  0.0, 12.0 },
            {  0.0, 44.0, 45.0, 34.0,  0.0 },
            {  0.0, 45.0, 55.0,  0.0,  0.0 },
            {  0.0, 34.0,  0.0, 33.0, 23.0 },
            { 12.0,  0.0,  0.0, 23.0, 24.0 }
        };

        public static readonly double[] cscValues = { 11.0, 12.0, 44.0, 45.0, 34.0, 45.0, 55.0, 34.0, 33.0, 23.0, 12.0, 23.0, 22.0 };
        public static readonly int[] cscRowIndices = {   0,    4,    1,    2,    3,    1,    2,    1,    3,    4,   0,     3,    4 };
        public static readonly int[] cscColOffsets = { 0, 2, 5, 7, 10, 13 };

        public static readonly double[] cscUpperValues = { 11.0, 44.0, 45.0, 55.0, 34.0, 33.0, 12.0, 23.0, 22.0 };
        public static readonly int[] cscUpperRowIndices = {   0,    1,    1,    2,    1,    3,    0,    3,    4 };
        public static readonly int[] cscUpperColOffsets = {   0, 1, 2, 4, 6, 9 };

        public static readonly int[] matlabPermutationAMD = { 2, 1, 3, 0, 4 }; // in 0-based indexing ofc
    }
}
