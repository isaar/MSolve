using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Commons
{
    //TODO: add precision, justification
    public class Array1DFormatting
    {
        public static readonly Array1DFormatting Default = new Array1DFormatting ();

        public Array1DFormatting(string separator = " ", string start="", string end="")
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
