using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationProvider : IAnalyzerProvider
    {
        void CalculateEffectiveMatrix(ISolverSubdomain subdomain, ImplicitIntegrationCoefficients coefficients);
        void ProcessRHS(ISolverSubdomain subdomain, ImplicitIntegrationCoefficients coefficients);
        void GetRHSFromHistoryLoad(int timeStep);
        void MassMatrixVectorProduct(ISolverSubdomain subdomain, IVector<double> vIn, double[] vOut);
        void DampingMatrixVectorProduct(ISolverSubdomain subdomain, IVector<double> vIn, double[] vOut);
    }
}
