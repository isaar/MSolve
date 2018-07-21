using System;
using System.Collections.Generic;
using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads 1D arrays from files.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array1DReader
    {
        private readonly bool firstLineIsLength;
        private readonly char[] separators;

        /// <summary>
        /// Initializes a new instance of the <see cref="Array1DReader"/> class, with the provided settings.
        /// </summary>
        /// <param name="firstLineIsLength">If true, the first line will be treated as the length of the array to be read. 
        ///     This allows optimizations, but will fail if there is no such line.</param>
        /// <param name="separators">The file is read as one string, which is then split whenever one of the characters in 
        ///     <paramref name="separators"/> is met. If none is provided, then the string will be split at each whitespace
        ///     character.</param>
        public Array1DReader(bool firstLineIsLength, params char[] separators)
        {
            this.firstLineIsLength = firstLineIsLength;
            this.separators = separators;
        }

        /// <summary>
        /// Reads the file at <paramref name="path"/> and returns an 1D array.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the length of the array specified in the first 
        ///     line is invalid.</exception>
        public double[] ReadFile(string path)
        {
            //TODO: add input checking
            if (firstLineIsLength) return ReadKnownSize(path);
            else return ReadUnknownSize(path);
        }

        private double[] ReadKnownSize(string path)
        {
            using (var reader = new StreamReader(path))
            {
                // Read the array length
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                int length = Int32.Parse(firstLine);
                var array = new double[length];

                // Read the array entries.
                string line = reader.ReadToEnd();
                string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                if (subStrings.Length != length) throw new IOException("Mismatch in provided entries and their declared count.");
                for (int i = 0; i < length; ++i)
                {
                    array[i] = Double.Parse(subStrings[i]);
                }

                return array;
            }
        }

        private double[] ReadUnknownSize(string path)
        {
            using (var reader = new StreamReader(path))
            {
                var array = new List<double>();
                string line = reader.ReadToEnd();
                string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < subStrings.Length; ++i)
                {
                    array.Add(Double.Parse(subStrings[i]));
                }
                return array.ToArray();
            }
        }
    }
}
