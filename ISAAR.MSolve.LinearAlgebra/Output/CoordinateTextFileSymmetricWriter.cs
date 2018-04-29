using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class CoordinateTextFileSymmetricWriter : MatrixWriter
    {
        private readonly ISparseSymmetricMatrix matrix;

        public CoordinateTextFileSymmetricWriter(ISparseSymmetricMatrix matrix)
        {
            this.matrix = matrix;
        }

        public static INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        protected override void WriteToStream(StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write($"{matrix.NumRows} {matrix.NumColumns} {matrix.CountNonZerosSuperDiagonal()}");
            foreach (var (row, col, val) in matrix.EnumerateNonZerosSuperDiagonal())
            {
                writer.WriteLine();
                writer.Write($"{row} {col} ");
                writer.Write(numberFormat, val);
            }
        }
    }
}
