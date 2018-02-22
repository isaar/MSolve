using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ILinearSystem
    {
        int ID { get; }
        IMatrix2D Matrix { get; set; }
        IVector RHS { get; }
        IVector Solution { get; set;  }
        //IVector GetRHSFromSolution(IVector solution, IVector dSolution);
        //void SubdomainToGlobalVector(double[] vIn, double[] vOut);
        //void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut);
        //void SplitGlobalVectorToSubdomain(double[] vIn, double[] vOut);
        //void SaveMaterialState();
        //void ClearMaterialStresses();
        //// REMOVE
        //void CloneMatrix();
    }
}
