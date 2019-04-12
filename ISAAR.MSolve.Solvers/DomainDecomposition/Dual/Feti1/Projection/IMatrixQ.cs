using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs a better name. What is the physical meaning of Q?
//TODO: Also implement the superlumped Q matrix.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Projection
{
    /// <summary>
    /// Represents the Q matrix of the FETI-1 projection operator.
    /// </summary>
    internal interface IMatrixQ
    {
        Vector Multiply(Vector vector);
        Matrix Multiply(Matrix matrix);
    }
}
