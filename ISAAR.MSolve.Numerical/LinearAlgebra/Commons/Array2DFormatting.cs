using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Commons
{
    //TODO: add precision, justification
    public class Array2DFormatting
    {
        public static readonly Array2DFormatting Default = new Array2DFormatting();

        public Array2DFormatting(string rowSeparator = "\n", string colSeparator = " ", string rowStart = "", string rowEnd = "")
        {
            this.RowSeparator = rowSeparator;
            this.ColSeparator = colSeparator;
            this.RowStart = rowStart;
            this.RowEnd = rowEnd;
        }

        public string ColSeparator { get; }
        public string RowSeparator { get; }
        public string RowStart { get; }
        public string RowEnd { get; }
    }
}
