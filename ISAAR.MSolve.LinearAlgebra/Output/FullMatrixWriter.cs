using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class FullMatrixWriter: MatrixWriter
    {
        private readonly Array2DFormatting arrayFormat;
        private readonly IIndexable2D matrix;

        public FullMatrixWriter(IIndexable2D matrix, Array2DFormatting format = null)
        {
            this.arrayFormat = (format == null) ? Array2DFormatting.Plain : format;
            this.matrix = matrix;
        }

        public static INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        protected override void WriteToStream(StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write(arrayFormat.ArrayStart);

            // First row
            writer.Write(arrayFormat.RowSeparator + arrayFormat.RowStart);
            writer.Write(string.Format(numberFormat, matrix[0, 0]));
            for (int j = 1; j < matrix.NumColumns; ++j)
            {
                writer.Write(arrayFormat.ColSeparator + string.Format(numberFormat, matrix[0, j]));
            }
            writer.Write(arrayFormat.RowEnd);

            // Subsequent rows
            for (int i = 1; i < matrix.NumRows; ++i)
            {
                writer.Write(arrayFormat.RowSeparator + arrayFormat.RowStart);
                writer.Write(string.Format(numberFormat, matrix[i, 0]));
                for (int j = 1; j < matrix.NumColumns; ++j)
                {
                    writer.Write(arrayFormat.ColSeparator + string.Format(numberFormat, matrix[i, j]));
                }
                writer.Write(arrayFormat.RowEnd);
            }
            writer.Write(arrayFormat.RowSeparator + arrayFormat.ArrayEnd);
        }
    }
}
