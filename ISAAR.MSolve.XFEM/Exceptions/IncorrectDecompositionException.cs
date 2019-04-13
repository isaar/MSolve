using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Exceptions
{
    class IncorrectDecompositionException: Exception
    {
        public IncorrectDecompositionException()
        { }

        public IncorrectDecompositionException(string message) : base(message)
        { }

        public IncorrectDecompositionException(string message, Exception inner) : base(message, inner)
        { }
    }
}
