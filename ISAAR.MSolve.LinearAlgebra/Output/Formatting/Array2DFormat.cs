using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    //TODO: add precision, justification
    public class Array2DFormat
    {
        public static readonly Array2DFormat Brackets = new Array2DFormat("\n[", "]\n", "[ ", " ]");
        public static readonly Array2DFormat Plain = new Array2DFormat("", "", "", "");

        public Array2DFormat(string arrayStart, string arrayEnd, string rowStart, string rowEnd, 
            string rowSeparator = "\n", string colSeparator = " ")
        {
            this.ArrayStart = arrayStart;
            this.ArrayEnd = arrayEnd;
            this.RowSeparator = rowSeparator;
            this.ColSeparator = colSeparator;
            this.RowStart = rowStart;
            this.RowEnd = rowEnd;
        }

        public string ArrayStart { get; }
        public string ArrayEnd { get; }
        public string ColSeparator { get; }
        public string RowSeparator { get; }
        public string RowStart { get; }
        public string RowEnd { get; }
    }
}
