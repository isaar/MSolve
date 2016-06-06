using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Solvers.PCGSkyline
{
    public class SolverPCGMatrixInitialization<T> : ISolverPCGInitialization where T : IMatrix2D<double>
    {
        private readonly SolverPCG<T> solver;

        public SolverPCGMatrixInitialization(SolverPCG<T> solver)
        {
            this.solver = solver;
        }

        public double InitializeAndGetResidual(IList<ISolverSubdomain> subdomains, double[] r, double[] x)
        {
            if (subdomains.Count != 1) throw new InvalidOperationException("Skyline PCG solver operates on one subdomain only.");

            double detf = 0;
            double temp = 0;
            ISolverSubdomain subdomain = subdomains[0];

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
