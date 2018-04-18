using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    class NonMatchingDimensionsException: ArgumentException
    {
        public NonMatchingDimensionsException(string message) : base(message) { }
    }
}
