using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class FullVectorWriter: MatrixWriter
    {
        private readonly Array1DFormatting format;
        private readonly bool writeLengthFirst = false;
        private readonly IVectorView vector;
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="firstLineIsLength">If true, the length of the vector will be written in the first line and then the 
        ///     vector entries will follow starting from the second.</param>
        /// <param name="format"></param>
        public FullVectorWriter(IVectorView vector, bool firstLineIsLength = false, Array1DFormatting format = null)
        {
            this.format = (format == null) ? Array1DFormatting.Plain : format;
            this.writeLengthFirst = firstLineIsLength;
            this.vector = vector;
        }

        public static INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        protected override void WriteToStream(StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            string separator = format.Separator;
            if (writeLengthFirst) writer.WriteLine(vector.Length);
            writer.Write(format.Start);
            writer.Write(string.Format(numberFormat, vector[0]));
            for (int i = 1; i < vector.Length; ++i)
            {
                writer.Write(separator + string.Format(numberFormat, vector[i]));
            }
            writer.WriteLine(format.End);
        }
    }
}
