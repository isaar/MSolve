using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning
{
    /// <summary>
    /// Implements the null object pattern in the contect of preconditioning. Use this class if you want to pass an 
    /// <see cref="IPreconditioner"/> object without actually applying any preconditioning, e.g. for benchmarking an iterative  
    /// algorithm. Using this preconditioner with PCG is equivalent to using CG, however the computational cost will be higher,
    /// since the operation z = inv(M) * r cannot be safely avoided; it just reduces to a vector copy.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class IdentityPreconditioner: IPreconditioner
    {
        /// <summary>
        /// Initializes a new instance of <see cref="IdentityPreconditioner"/> with the provided settings.
        /// </summary>
        public IdentityPreconditioner()
        {
        }

        /// <summary>
        /// See <see cref="IPreconditioner.SolveLinearSystem(Vector)"/>.
        /// </summary>
        /// <remarks>
        /// This method works for all dimensions of the preconditioner matrix and the right hand side vector. This way the user
        /// doesn't have to define the dimensions of the linear system, which is useful when testing or benchmarking, at the 
        /// expense of little extra safety.
        /// </remarks>
        public void SolveLinearSystem(IVectorView rhsVector, IVector lhsVector) => lhsVector.CopyFrom(rhsVector);

        /// <summary>
        /// Creates instances of <see cref="IdentityPreconditioner"/>.
        /// </summary>
        public class Factory: IPreconditionerFactory
        {
            /// <summary>
            /// Initializes a new instance of <see cref="IdentityPreconditioner.Factory"/>.
            /// </summary>
            public Factory() { }

            /// <summary>
            /// See <see cref="IPreconditionerFactory.CreatePreconditionerFor(IMatrixView)"/>.
            /// </summary>
            public IPreconditioner CreatePreconditionerFor(IMatrixView matrix) => new IdentityPreconditioner();
        }
    }
}
