using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

// Inversion is best handled by the matrix object itself, since the original should overwrite the factorized data in most cases, 
// which should be hidden from the user. Besides, I am not sure if first factorizing the matrix is more efficient than 
// Gauss-Jordan.
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    public interface IFactorization
    {
        double CalcDeterminant();
        Vector SolveLinearSystem(Vector rhs);
    }
}
