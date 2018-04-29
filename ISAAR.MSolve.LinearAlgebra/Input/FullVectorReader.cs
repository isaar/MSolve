using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    public class FullVectorReader
    {
        private static readonly char[] separators = { };
        private readonly bool firstLineIsLength;

        public FullVectorReader(bool firstLineIsLength)
        {
            this.firstLineIsLength = firstLineIsLength;
        }

        public Vector ReadFromFile(string path)
        {
            //TODO: add input checking
            if (firstLineIsLength) return Vector.CreateFromArray(ReadKnownSize(path), false);
            else return Vector.CreateFromArray(ReadUnknownSize(path), false);
        }

        private double[] ReadKnownSize(string path)
        {
            using (var reader = new StreamReader(path))
            {
                // Read the vector length
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                int length = Int32.Parse(firstLine);
                var vector = new double[length];

                // Read the vector entries. TODO: what if each entry is on a different line?
                string line = reader.ReadLine();
                string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                if (subStrings.Length != length) throw new IOException("Mismatch in provided entries and their declared count.");
                for (int i = 0; i < length; ++i)
                {
                    vector[i] = Double.Parse(subStrings[i]);
                }

                return vector;
            }
        }

        private double[] ReadUnknownSize(string path)
        {
            using (var reader = new StreamReader(path))
            {
                var vector = new List<double>();
                string line = reader.ReadLine(); //TODO: what if each entry is on a different line?
                string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < subStrings.Length; ++i)
                {
                    vector.Add(Double.Parse(subStrings[i]));
                }
                return vector.ToArray();
            }
        }
    }
}
