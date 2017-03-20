using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class ValueComparer
    {
        private readonly double tolerance;

        public ValueComparer(double tolerance)
        {
            if (tolerance <= 0.0 || tolerance > 1.0)
                throw new ArgumentException("tolerance must belong to the interval (0,1]");
            this.tolerance = tolerance;
        }

        public bool AreEqual(double actual, double expected)
        {
            if (Math.Abs(expected) < tolerance && Math.Abs(actual) < tolerance) return true;
            else if (Math.Abs(expected) < tolerance && Math.Abs(actual) >= tolerance) return false;
            else return (Math.Abs(1.0 - actual / expected) < tolerance) ? true : false;
        }
    }
}
