using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.PSM
{
    public class SolverPSMMatrixInitialization<T> : ISolverPCGInitialization where T : IMatrix2D
    {
        private readonly SolverPSM solver;

        public SolverPSMMatrixInitialization(SolverPSM solver)
        {
            this.solver = solver;
        }

        public double InitializeAndGetResidual(IList<ILinearSystem> subdomains, IVector r, IVector x)
        {
            if (subdomains.Count <= 1) throw new InvalidOperationException("PSM solver operates on more than one subdomains.");

            
            double detf = 0;
            double temp = 0;
            ILinearSystem subdomain = subdomains[0];

            for (int i = 0; i < subdomain.RHS.Length; i++)
            {
                temp = subdomain.RHS[i];
                detf += temp * temp;
                r[i] = temp;
            }

            return Math.Sqrt(detf);
        }
    }
}
