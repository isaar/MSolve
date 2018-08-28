using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning
{
    /// <summary>
    /// Implements the null object pattern in the contect of preconditioning. Use this class if you want to pass an 
    /// <see cref="IPreconditioner"/> object without actually applying any preconditioning, e.g. for benchmarking an iterative  
    /// algorithm. It can avoid many performance costs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class IdentityPreconditioner: IPreconditioner
    {
        private readonly bool copyRhs;

        /// <summary>
        /// Initializes a new instance of <see cref="IdentityPreconditioner"/> with the provided settings.
        /// </summary>
        /// <param name="copyRhs">True to copy the rhs vector during <see cref="SolveLinearSystem(Vector)"/> and ensure it is 
        ///     not overwritten by the iterative algorithm. False to avoid the overhead of copying the rhs vector.</param>
        public IdentityPreconditioner(bool copyRhs = true)
        {
            this.copyRhs = copyRhs;
        }

        /// <summary>
        /// See <see cref="IPreconditioner.SolveLinearSystem(Vector)"/>.
        /// </summary>
        /// <remarks>
        /// This method works for all dimensions of the preconditioner matrix and the right hand side vector. This way the user
        /// doesn't have to define the dimensions of the linear system, which is useful when testing or benchmarking, at the 
        /// expense of little extra safety.
        /// </remarks>
        public Vector SolveLinearSystem(Vector rhs)
        {
            if (copyRhs) return rhs.Copy();
            else return rhs;
        }
    }
}
