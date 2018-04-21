using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Preconditioning
{
    public interface IPreconditioner
    {
        /// <summary>
        /// Apply the preconditioner. Solves the system M * v = w, where M is the preconditioner and the definition of the vectors 
        /// v, w depends on the iterative algorithm. 
        /// </summary>
        /// <param name="rhs">The right hand side vector of the system M * v = w. It may or may not be the right hand side vector
        ///     of the original linear system.</param>
        /// <returns></returns>
        Vector SolveLinearSystem(Vector rhs);
    }
}
