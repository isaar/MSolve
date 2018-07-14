using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Rename to PatternModified or something else that also works for vectors.
namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    public class MatrixPatternException: Exception
    {
        public MatrixPatternException()
        { }

        public MatrixPatternException(string message) : base(message)
        { }

        public MatrixPatternException(string message, Exception inner) : base(message, inner)
        { }
    }
}
