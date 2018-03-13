using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Optimization.Initialization.HaltonPoints
{
    interface IPrimeGenerator
    {
        int[] FirstPrimes(int count);
    }
}
