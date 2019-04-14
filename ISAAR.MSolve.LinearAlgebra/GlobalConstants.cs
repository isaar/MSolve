using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra
{
    /// <summary>
    /// Constant values that are not provided by the standard library.
    /// Author: Serafeim Bakalakos
    /// </summary>
    public static class GlobalConstants
    {
        static GlobalConstants()
        {
            MachinePrecisionDouble = 1.0d;
            do MachinePrecisionDouble /= 2.0d;
            while ((double)(1.0 + MachinePrecisionDouble) != 1.0);
        }

        public static readonly double MachinePrecisionDouble;
    }
}
