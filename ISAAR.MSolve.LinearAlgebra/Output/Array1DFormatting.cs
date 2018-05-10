using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    //TODO: add precision, justification
    public class Array1DFormatting
    {
        public static readonly Array1DFormatting Brackets = new Array1DFormatting("[ ", " ]");
        public static readonly Array1DFormatting PlainHorizontal = new Array1DFormatting("", "", " ");
        public static readonly Array1DFormatting PlainVertical = new Array1DFormatting("", "", "\n");

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
