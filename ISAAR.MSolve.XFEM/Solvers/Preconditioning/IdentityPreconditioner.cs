using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Solvers.Preconditioning
{
    /// <summary>
    /// Use this class if you do not want to apply any preconditioning. It works for all matrix and vector dimensions. It can
    /// avoid many performance costs.
    /// </summary>
    public class IdentityPreconditioner: IPreconditioner
    {
        private readonly bool copyRhs;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="copyRhs">True to copy the rhs vector and ensure it is not overwritten by the iterative algorithm. 
        ///     False to avoid the overhead of copying the rhs vector.</param>
        public IdentityPreconditioner(bool copyRhs = true)
        {
            this.copyRhs = copyRhs;
        }

        public Vector SolveLinearSystem(Vector rhs)
        {
            if (copyRhs) return rhs.Copy();
            else return rhs;
        }
    }
}
