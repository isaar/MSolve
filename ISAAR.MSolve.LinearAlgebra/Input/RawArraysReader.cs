using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads the internal indexing and value arrays of a sparse matrix from a file.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RawArraysReader
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="RawArraysReader"/> class.
        /// </summary>
        public RawArraysReader()
        { }

        /// <summary>
        /// If false, an exception will be thrown if an empty line is met during reading the file.
        /// </summary>
        private bool IgnoreEmptyLines { get; set; } = true;

        /// <summary>
        /// Reads the values and diagonal offsets arrays of a matrix in Skyline format from 2 separate files.
        /// </summary>
        /// <param name="valuesArrayPath">Each entry on a separate line. The first must be the array length.</param>
        /// <param name="diagOffsetsArrayPath">Each entry on a separate line. The first must be the array length.</param>
        /// <param name="checkSkylineArrays">
        /// If true, the provided arrays will be checked to make sure they are valid Skyline arrays, which is safer. 
        /// If false, no such check will take place, which is faster.
        /// </param>
        public SkylineMatrix ReadSkylineMatrixFromSeparateFiles(string valuesArrayPath, string diagOffsetsArrayPath, 
            bool checkSkylineArrays)
        {
            // Values array
            double[] values;
            using (var reader = new StreamReader(valuesArrayPath))
            {
                string line = reader.ReadLine();
                int length = int.Parse(line.Trim());
                values = new double[length];
                for (int i = 0; i < length; ++i)
                {
                    line = reader.ReadLine();
                    values[i] = double.Parse(line.TrimEnd());
                }
            }

            // Diagonal offsets array
            int[] diagOffsets;
            using (var reader = new StreamReader(diagOffsetsArrayPath))
            {
                string line = reader.ReadLine();
                int length = int.Parse(line.Trim());
                diagOffsets = new int[length];
                for (int i = 0; i < length; ++i)
                {
                    line = reader.ReadLine();
                    diagOffsets[i] = int.Parse(line.TrimEnd());
                }
            }

            return SkylineMatrix.CreateFromArrays(diagOffsets.Length - 1, values, diagOffsets, checkSkylineArrays);
        }
    }
}
