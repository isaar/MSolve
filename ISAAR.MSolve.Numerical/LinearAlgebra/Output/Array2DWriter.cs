using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    class Array2DWriter: MatrixWriter
    {
        private readonly Array2DFormatting format;
        private readonly double[,] matrix;

        public Array2DWriter(double[,] matrix, Array2DFormatting format = null)
        {
            this.format = (format == null) ? Array2DFormatting.Plain : format;
            this.matrix = matrix;
        }

        public static INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        protected override void WriteToStream(StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write(format.ArrayStart);

            // First row
            writer.Write(format.RowSeparator + format.RowStart);
            writer.Write(string.Format(numberFormat, matrix[0, 0]));
            for (int j = 1; j < matrix.GetLength(1); ++j)
            {
                writer.Write(format.ColSeparator + string.Format(numberFormat, matrix[0, j]));
            }
            writer.Write(format.RowEnd);

            // Subsequent rows
            for (int i = 1; i < matrix.GetLength(0); ++i)
            {
                writer.Write(format.RowSeparator + format.RowStart);
                writer.Write(string.Format(numberFormat, matrix[i, 0]));
                for (int j = 1; j < matrix.GetLength(1); ++j)
                {
                    writer.Write(format.ColSeparator + string.Format(numberFormat, matrix[i, j]));
                }
                writer.Write(format.RowEnd);
            }
            writer.Write(format.RowSeparator + format.ArrayEnd);
        }
    }
}
