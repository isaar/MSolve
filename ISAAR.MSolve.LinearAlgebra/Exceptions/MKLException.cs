using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    public class MKLException : Exception
    {
        public MKLException()
        { }

        public MKLException(string message) : base(message)
        { }

        public MKLException(string message, Exception inner) : base(message, inner)
        { }
    }
}
