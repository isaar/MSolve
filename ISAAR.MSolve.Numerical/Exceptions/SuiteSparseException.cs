using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.Exceptions
{
    public class SuiteSparseException : Exception
    {
        public SuiteSparseException()
        { }

        public SuiteSparseException(string message) : base(message)
        { }

        public SuiteSparseException(string message, Exception inner) : base(message, inner)
        { }
    }
}
