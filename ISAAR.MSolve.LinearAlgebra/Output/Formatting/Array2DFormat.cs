using System;

namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    /// <summary>
    /// Describes how entries of a 2D array or matrix are separated during output.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array2DFormat
    {
        /// <summary>
        /// E.g.
        /// [
        /// [ 1 2 3 ]
        /// [ 4 5 6 ]
        /// ]
        /// 
        /// </summary>
        public static readonly Array2DFormat Brackets = new Array2DFormat("\n[", "]\n", "[ ", " ]" + Environment.NewLine, " ");

        /// <summary>
        /// E.g. (notice the comma after the last row)
        /// {
        ///     { 1 2 3 },
        ///     { 4 5 6 },
        /// }
        /// </summary>
        public static readonly Array2DFormat CSharpArray = new Array2DFormat("{", "}" + Environment.NewLine, "\t{ ",
            " }," + Environment.NewLine, ", ");

        /// <summary>
        /// E.g.
        /// 1 2 3
        /// 4 5 6
        /// </summary>
        public static readonly Array2DFormat Plain = new Array2DFormat("", "", "", Environment.NewLine, " ");
        

        /// <summary>
        /// Initializes a new instance of the <see cref="Array2DFormat"/> with the provided properties.
        /// </summary>
        /// <param name="arrayStart">A string to write before everything else.</param>
        /// <param name="arrayEnd">A string to write after everything else.</param>
        /// <param name="rowStart">A string to write before the entries of each row.</param>
        /// <param name="rowSeparator">A string to write between two consecutive rows.</param>
        /// <param name="colSeparator">A string to write between two consecutive entries of the same row.</param>
        public Array2DFormat(string arrayStart, string arrayEnd, string rowStart, string rowSeparator, string colSeparator = " ") //TODO: rowSeparator and rowEnd should be one
        {
            this.ArrayStart = arrayStart;
            this.ArrayEnd = arrayEnd;
            this.RowSeparator = rowSeparator;
            this.ColSeparator = colSeparator;
            this.RowStart = rowStart;
        }

        /// <summary>
        /// A string to write after everything else.
        /// </summary>
        public string ArrayEnd { get; }

        /// <summary>
        /// A string to write before everything else.
        /// </summary>
        public string ArrayStart { get; }

        /// <summary>
        /// A string to write between two consecutive entries of the same row.
        /// </summary>
        public string ColSeparator { get; }

        /// <summary>
        /// A string to write after the entries of each row.
        /// </summary>
        public string RowEnd { get; }

        /// <summary>
        /// A string to write between two consecutive rows.
        /// </summary>
        public string RowSeparator { get; }

        /// <summary>
        /// A string to write before the entries of each row.
        /// </summary>
        public string RowStart { get; }
    }
}
