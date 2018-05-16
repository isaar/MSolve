using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.LinearAlgebra.Exceptions
{
    class SymmetricPatternModifiedException : MatrixPatternException
    {
        public SymmetricPatternModifiedException()
        { }

        public SymmetricPatternModifiedException(string message) : base(message)
        { }

        public SymmetricPatternModifiedException(string message, Exception inner) : base(message, inner)
        { }
    }
}
