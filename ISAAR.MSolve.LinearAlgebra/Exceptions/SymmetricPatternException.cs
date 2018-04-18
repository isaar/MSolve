using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    class SymmetricPatternException : MatrixPatternException
    {
        public SymmetricPatternException()
        { }

        public SymmetricPatternException(string message) : base(message)
        { }

        public SymmetricPatternException(string message, Exception inner) : base(message, inner)
        { }
    }
}
