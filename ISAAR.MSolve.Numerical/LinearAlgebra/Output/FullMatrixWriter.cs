using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    public class FullMatrixWriter: MatrixWriter
    {
        private readonly Array2DFormatting format;
        private readonly IIndexable2D matrix;

        public FullMatrixWriter(IIndexable2D matrix, Array2DFormatting format = null)
        {
            this.format = (format == null) ? Array2DFormatting.Plain : format;
            this.matrix = matrix;
        }

        protected override void WriteToStream(StreamWriter writer)
        {
            writer.Write(format.ArrayStart);

            // First row
            writer.Write(format.RowSeparator + format.RowStart);
            writer.Write(matrix[0, 0]);
            for (int j = 1; j < matrix.NumColumns; ++j)
            {
                writer.Write(format.ColSeparator + matrix[0, j]);
            }
            writer.Write(format.RowEnd);

            // Subsequent rows
            for (int i = 1; i < matrix.NumRows; ++i)
            {
                writer.Write(format.RowSeparator + format.RowStart);
                writer.Write(matrix[i, 0]);
                for (int j = 1; j < matrix.NumColumns; ++j)
                {
                    writer.Write(format.ColSeparator + matrix[i, j]);
                }
                writer.Write(format.RowEnd);
            }
            writer.Write(format.RowSeparator + format.ArrayEnd);
        }
    }
}
