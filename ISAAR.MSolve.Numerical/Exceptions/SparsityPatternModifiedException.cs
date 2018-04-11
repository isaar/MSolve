using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.Exceptions
{
    public class SparsityPatternModifiedException: MatrixPatternException
    {
        public SparsityPatternModifiedException()
        { }

        public SparsityPatternModifiedException(string message) : base(message)
        { }

        public SparsityPatternModifiedException(string message, Exception inner) : base(message, inner)
        { }
    }
}
