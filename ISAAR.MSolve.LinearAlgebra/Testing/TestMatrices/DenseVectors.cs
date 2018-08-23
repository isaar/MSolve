using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class DenseVectors
    {
        public const int length = 10;

        public static readonly double[] vector1 = {
            4.602954, 4.717625, 0.344063, 3.068357, 1.371406, 4.356552, 1.298073, 1.672528, 2.855245, 0.056439 };

        public static void Print()
        {
            var vector = Vector.CreateFromArray(vector1);
            Console.WriteLine("Vector = ");
            var writer = new FullVectorWriter(false) { ArrayFormat = Array1DFormat.Brackets };
            writer.WriteToConsole(vector);
        }
    }
}
