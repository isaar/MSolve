namespace ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces
{
    public interface ISolveable
    {
        void Solve(IVectorOLD f, IVectorOLD result);
    }
}
