using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.Commons
{
    public class SingleSubdomainSystem<TMatrix> : LinearSystem_v2<TMatrix, Vector>
        where TMatrix : class, IMatrixView //TODO: perhaps this should be IMatrix
    {
        internal SingleSubdomainSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
        internal override Vector CreateZeroVector() => Vector.CreateZero(this.Size);
    }
}
