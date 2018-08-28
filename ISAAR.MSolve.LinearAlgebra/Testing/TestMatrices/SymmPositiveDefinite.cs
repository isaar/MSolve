using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    public class SymmPositiveDefinite
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            { 8.9156,    0.4590,    0.0588,    0.5776,    0.7118,    0.7423,    0.4389,    0.4353,    0.4929,    0.7223 },
            { 0.4590,    7.5366,    0.4276,    0.3282,    0.5277,    0.4274,    0.2498,    0.7622,    0.4987,    0.8953 },
            { 0.0588,    0.4276,    6.4145,    0.4144,    0.5954,    0.6196,    0.3257,    0.5084,    0.6342,    0.5270 },
            { 0.5776,    0.3282,    0.4144,    7.7576,    0.6609,    0.9212,    0.8040,    0.2146,    0.5077,    0.3928 },
            { 0.7118,    0.5277,    0.5954,    0.6609,    9.0882,    0.5096,    0.3434,    0.6139,    0.7590,    0.5276 },
            { 0.7423,    0.4274,    0.6196,    0.9212,    0.5096,   10.5523,    0.4885,    0.7861,    0.5294,    0.7933 },
            { 0.4389,    0.2498,    0.3257,    0.8040,    0.3434,    0.4885,    7.5277,    0.6574,    0.4020,    0.2763 },
            { 0.4353,    0.7622,    0.5084,    0.2146,    0.6139,    0.7861,    0.6574,    7.5142,    0.3528,    0.6254 },
            { 0.4929,    0.4987,    0.6342,    0.5077,    0.7590,    0.5294,    0.4020,    0.3528,    7.8487,    0.3035 },
            { 0.7223,    0.8953,    0.5270,    0.3928,    0.5276,    0.7933,    0.2763,    0.6254,    0.3035,    9.6122 }};

        public static readonly double[] lhs = { 2.6621, 3.5825, 0.8965, 1.6827, 0.9386, 1.6096, 2.0193, 2.7428, 0.2437, 2.7637 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static void Print()
        {
            var A = SymmetricMatrix.CreateFromArray(matrix);
            Console.WriteLine("Symmetric positive definite matrix = ");
            var writer = new FullMatrixWriter(); // TODO: implement triangular printer
            writer.WriteToConsole(A);
        }
    }
}
