using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    //TODO: add precision, justification
    public class Array1DFormatting
    {
        public static readonly Array1DFormatting Brackets = new Array1DFormatting("[ ", " ]");
        public static readonly Array1DFormatting Plain = new Array1DFormatting("", "");

        public Array1DFormatting(string start, string end, string separator = " ")
        {
            this.Separator = separator;
            this.Start = start;
            this.End = end;
        }

        public string End { get; }
        public string Separator { get; }
        public string Start { get; }
    }
}
