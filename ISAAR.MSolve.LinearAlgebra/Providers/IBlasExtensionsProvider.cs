using System;
using System.Collections.Generic;
using System.Text;

//TODO: this should be internal
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public interface IBlasExtensionsProvider
    {
        /// <summary>
        /// y = alpha * x + beta * y
        /// </summary>
        void Daxpby(int n, double alpha, double[] x, int offsetX, int incX, double beta, double[] y, int offsetY, int incY);
    }
}
