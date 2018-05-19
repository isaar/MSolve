using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    class EmptyMatrixBuilderException: Exception
    {
        public EmptyMatrixBuilderException()
        { }

        public EmptyMatrixBuilderException(string message) : base(message)
        { }

        public EmptyMatrixBuilderException(string message, Exception inner) : base(message, inner)
        { }
    }
}
