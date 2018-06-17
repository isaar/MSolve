using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class ElasticityTensorContinuum2D : Matrix2D
    {
        public ElasticityTensorContinuum2D() : base(3, 3) { }

        public ElasticityTensorContinuum2D(double[,] data) : base(data)
        {
            if (data.GetLength(0) != 3 || data.GetLength(1) != 3)
                throw new ArgumentException($"input array (dimensions: {data.GetLength(0)} by {data.GetLength(1)}) is not 3 by 3");
        }
    }
}
