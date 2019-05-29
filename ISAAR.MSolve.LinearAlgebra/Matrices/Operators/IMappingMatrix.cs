using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices.Operators
{
    public interface IMappingMatrix : IIndexable2D
    {
        Vector Multiply(Vector vector, bool transposeThis = true);
    }
}
