using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.Interfaces
{
    public interface ISearchVectorCalculator
    {
        bool InitializeStartingVectorFromSearchVectors(IVector<double> x, IVector<double> b);
        void CalculateSearchVector(IIterativeSolver solver);
        double CalculateGradient(IIterativeSolver solver);
        void ClearSearchVectors(int vectorsToKeepFromTop);
    }
}
