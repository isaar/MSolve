namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface ISolveable
    {
        void Solve(IVector f, IVector result);
    }
}
