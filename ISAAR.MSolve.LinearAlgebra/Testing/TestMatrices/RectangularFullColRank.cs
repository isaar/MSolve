using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class RectangularFullColRank
    {
        public const int numRows = 10;
        public const int numCols = 5;

        public static readonly double[,] matrix = new double[,] {
            { 1.687338, 1.910625, 0.078085, 2.332444, 3.131424 },
            { 3.215230, 1.767679, 1.883659, 0.809185, 2.026953 },
            { 0.981285, 1.325970, 2.639193, 0.923645, 3.475809 },
            { 0.751857, 2.296387, 2.101308, 1.871156, 3.148936 },
            { 3.446088, 3.932793, 0.828995, 0.343996, 1.349877 },
            { 0.779645, 1.764859, 1.694237, 2.854230, 1.260338 },
            { 0.968975, 3.746153, 2.681011, 2.329230, 2.481479 },
            { 3.581092, 1.647008, 3.430712, 2.409152, 3.780291 },
            { 1.536053, 0.897209, 3.872541, 1.485733, 0.601071 },
            { 3.186167, 0.066731, 1.536828, 2.245892, 2.542407 }};

        public static readonly double[] lhs5 = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064 };
        public static readonly double[] rhs5 = MatrixOperations.MatrixTimesVector(matrix, lhs5);
       
        public static void Print()
        {
            var A = Matrix.CreateFromArray(matrix);
            Console.WriteLine("Rectangular matrix = ");
            var writer = new FullMatrixWriter();
            writer.WriteToConsole(A);
        }
    }
}
