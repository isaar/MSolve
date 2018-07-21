using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Output.Formatting
{
    //TODO: add precision, justification
    public class Array1DFormat
    {
        public static readonly Array1DFormat Brackets = new Array1DFormat("[ ", " ]");
        public static readonly Array1DFormat PlainHorizontal = new Array1DFormat("", "", " ");
        public static readonly Array1DFormat PlainVertical = new Array1DFormat("", "", "\n");

        public Array1DFormat(string start, string end, string separator = " ")
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
