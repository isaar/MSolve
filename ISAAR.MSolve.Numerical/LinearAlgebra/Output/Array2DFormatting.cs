using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    //TODO: add precision, justification
    public class Array2DFormatting
    {
        public static readonly Array2DFormatting Brackets = new Array2DFormatting("\n[", "]\n", "[ ", " ]");
        public static readonly Array2DFormatting Plain = new Array2DFormatting("", "", "", "");

        public Array2DFormatting(string arrayStart, string arrayEnd, string rowStart, string rowEnd, 
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
