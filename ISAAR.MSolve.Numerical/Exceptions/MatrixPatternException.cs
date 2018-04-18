using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.Exceptions
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
