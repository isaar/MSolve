using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class ElasticityTensorContinuum3D : Matrix2D
    {
        public ElasticityTensorContinuum3D() : base(6, 6) { }

        public ElasticityTensorContinuum3D(double[,] data) : base(data)
        {
            if (data.GetLength(0) != 6 || data.GetLength(1) != 6)
                throw new ArgumentException($"input array (dimensions: {data.GetLength(0)} by {data.GetLength(1)}) is not 6 by 6");
        }
    }
}
