using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.PCGSkyline;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloSolverPCGDirectMatrixCalculator : ISolverPCGMatrixCalculator, ISolverPCGInitialization
    {
        private SolverPCG<SkylineMatrix2D> solver;
        private MonteCarloAnalyzerWithStochasticMaterial analyzer;
        private SkylineMatrix2D preconditioner;

        public SolverPCG<SkylineMatrix2D> Solver
        {
            get { return solver; }
            set { solver = value; }
        }

        public MonteCarloAnalyzerWithStochasticMaterial Analyzer
        {
            get { return analyzer; }
            set { analyzer = value; }
        }

        public int VectorSize
        {
            get
            {
                //return solver.SubdomainsDictionary.Values.First().RHS.Length;
                return solver.LinearSystem.RHS.Length;
            }
        }

        public void Precondition(IVector vIn, IVector vOut)
        {
            //if (analyzer.FactorizedMatrices.Count != 1)
            //    throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            //foreach (var m in analyzer.FactorizedMatrices.Values)
            //    m.Solve(vIn, ((Vector<double>)vOut).Data);
            preconditioner.Solve(vIn, vOut);
        }

        public void MultiplyWithMatrix(IVector vIn, IVector vOut)
        {
            //((SkylineMatrix2D)solver.SubdomainsDictionary.Values.First().Matrix).Multiply(vIn, ((Vector)vOut).Data, 1.0, 0, 0, true);
            ((SkylineMatrix2D)solver.LinearSystem.Matrix).Multiply(vIn, ((Vector)vOut).Data, 1.0, 0, 0, true);
        }

        public double InitializeAndGetResidual(IList<ILinearSystem> subdomains, IVector r, IVector x)
        {
            double detf = 0;
            double temp = 0;

            if (subdomains.Count != 1)
                throw new InvalidOperationException("Cannot initialize and calculate residuals with more than one subdomains");

            foreach (ILinearSystem subdomain in subdomains)
            {
                r.CopyFrom(0, subdomain.RHS.Length, subdomain.RHS, 0);
                //Array.Copy(((Vector)subdomain.RHS).Data, r, subdomain.RHS.Length);

                //subdomain.SubdomainToGlobalVector(((Vector<double>)subdomain.RHS).Data, r);
                var s = (SkylineMatrix2D)subdomain.Matrix;

                if (preconditioner == null)
                {
                    preconditioner = (SkylineMatrix2D)s.Clone();
                    preconditioner.Factorize(1e-8, new List<IVector>(), new List<int>());
                }
            }

            for (int i = 0; i < r.Length; i++)
            {
                temp = r[i];
                detf += temp * temp;
            }
            return Math.Sqrt(detf);
        }
    }
}
