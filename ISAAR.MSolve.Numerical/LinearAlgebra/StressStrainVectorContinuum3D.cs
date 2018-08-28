using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class StressStrainVectorContinuum3D : Vector
    {
        public StressStrainVectorContinuum3D() : base(6) { }
        public StressStrainVectorContinuum3D(double[] data) : base(data)
        {
            if (data.Length != 6)
                throw new ArgumentException($"input array dimension ({data.Length}) is not 6");
        }
    }
}
