using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class StressStrainVectorContinuum2D : Vector
    {
        public StressStrainVectorContinuum2D() : base(6) { }
        public StressStrainVectorContinuum2D(double[] data) : base(data)
        {
            if (data.Length != 3)
                throw new ArgumentException($"input array dimension ({data.Length}) is not 3");
        }
    }
}
