using System;

namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Describes how entries of a 1D array, list or vector are separated during output.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array1DFormat
    {
        /// <summary>
        /// E.g. 
        /// [ 1 2 3 ]
        /// </summary>
        public static readonly Array1DFormat Brackets = new Array1DFormat("[ ", " ]");

        /// <summary>
        /// E.g. 
        /// [ 1 2 3 ]
        /// </summary>
        public static readonly Array1DFormat CSharpArray = new Array1DFormat("{ ", " }", ", ");

        /// <summary>
        /// E.g. 
        /// 1 2 3
        /// </summary>
        public static readonly Array1DFormat PlainHorizontal = new Array1DFormat("", "", " ");

        /// <summary>
        /// E.g. 
        /// 1
        /// 2
        /// 3
        /// </summary>
        public static readonly Array1DFormat PlainVertical = new Array1DFormat("", "", Environment.NewLine);

        /// <summary>
        /// Initializes a new instance of the <see cref="Array1DFormat"/> with the provided properties.
        /// </summary>
        /// <param name="start">A string to write before all entries.</param>
        /// <param name="end">A string to write after all entries.</param>
        /// <param name="separator">A string to write between two consecutive entries.</param>
        public Array1DFormat(string start, string end, string separator = " ")
        {
            this.Separator = separator;
            this.Start = start;
            this.End = end;
        }

        /// <summary>
        /// A string to write after all entries.
        /// </summary>
        public string End { get; }

        /// <summary>
        /// A string to write between two consecutive entries.
        /// </summary>
        public string Separator { get; }

        /// <summary>
        /// A string to write written before all entries.
        /// </summary>
        public string Start { get; }
    }
}
