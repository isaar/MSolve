using System;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// Implementation of <see cref="ILinearSystem"/> that can be used with solvers withour domain decomposition.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix"></typeparam>
    public class SingleSubdomainSystem<TMatrix> : LinearSystemBase<TMatrix, Vector>, ISingleSubdomainLinearSystem
        where TMatrix : class, IMatrix
    {
        internal SingleSubdomainSystem(ISubdomain subdomain) : base(subdomain) { }

        public override IVector CreateZeroVector() => CreateZeroVectorConcrete();

        public Vector CreateZeroVectorConcrete()
        {
            if (Size == initialSize) throw new InvalidOperationException(
                "The linear system size must be set before creating vectors. First of all order the subdomain freedom degrees.");
            return Vector.CreateZero(Size);
        }        
    }
}
