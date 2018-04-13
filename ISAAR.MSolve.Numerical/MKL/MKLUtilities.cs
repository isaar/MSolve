using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.MKL
{
    // TODO: add a method that checks negative infos. It is basically the same everywhere in MKL. Only positive infos change 
    // depending on the function.
    static class MKLUtilities
    {
        /// <summary>
        /// Use this value to initialize the output arguments named info. So far, negative values of the info are parameter 
        /// indices, thus the int.MinValue is safe to use as a default.  
        /// </summary>
        public const int DefaultInfo = int.MinValue;
    }
}
