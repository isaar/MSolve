using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    public class Array1DWriter: MatrixWriter
    {
        private readonly Array1DFormatting format;
        private readonly double[] vector;

        public Array1DWriter(double[] vector, Array1DFormatting format = null)
        {
            this.format = (format == null) ? Array1DFormatting.Plain : format;
            this.vector = vector;
        }

        protected override void WriteToStream(StreamWriter writer)
        {
            string separator = format.Separator;
            writer.Write(format.Start);
            Console.Write(vector[0]);
            for (int i = 1; i < vector.Length; ++i)
            {
                writer.Write(separator + vector[i]);
            }
            writer.WriteLine(format.End);
        }
    }
}
