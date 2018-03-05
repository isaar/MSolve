using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.Exceptions
{
    class NonMatchingDimensionsException: Exception
    {
        public NonMatchingDimensionsException()
        { }

        public NonMatchingDimensionsException(string message): base(message)
        { }

        public NonMatchingDimensionsException(string message, Exception inner): base(message, inner)
        { }
    }
}
