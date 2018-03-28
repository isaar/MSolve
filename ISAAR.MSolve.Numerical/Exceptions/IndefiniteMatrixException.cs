using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.Exceptions
{
    class IndefiniteMatrixException : Exception
    {
        public IndefiniteMatrixException()
        { }

        public IndefiniteMatrixException(string message) : base(message)
        { }

        public IndefiniteMatrixException(string message, Exception inner) : base(message, inner)
        { }
    }
}
