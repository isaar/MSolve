using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    public class SingularMatrixException: IndefiniteMatrixException
    {
        public SingularMatrixException()
        { }

        public SingularMatrixException(string message): base(message)
        { }

        public SingularMatrixException(string message, Exception inner): base(message, inner)
        { }
    }
}
