using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISolverSubdomain
    {
        int ID { get; }
        IMatrix2D<double> Matrix { get; set; }
        IVector<double> RHS { get; }
        IVector<double> Solution { get; set;  }
        IVector<double> GetRHSFromSolution(IVector<double> solution, IVector<double> dSolution);
        void SubdomainToGlobalVector(double[] vIn, double[] vOut);
        void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut);
        void SplitGlobalVectorToSubdomain(double[] vIn, double[] vOut);
        void SaveMaterialState();
        void ClearMaterialStresses();
        //// REMOVE
        //void CloneMatrix();
    }
}
