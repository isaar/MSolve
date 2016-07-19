using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.PCGSkyline;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloSolverPCGDirectMatrixCalculator : ISolverPCGMatrixCalculator, ISolverPCGInitialization
    {
        private SolverPCG<SkylineMatrix2D<double>> solver;
        private MonteCarloAnalyzerWithStochasticMaterial analyzer;
        private SkylineMatrix2D<double> preconditioner;

        public SolverPCG<SkylineMatrix2D<double>> Solver
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
                return solver.SubdomainsDictionary.Values.First().RHS.Length;
            }
        }

        public void Precondition(IVector<double> vIn, IVector<double> vOut)
        {
            //if (analyzer.FactorizedMatrices.Count != 1)
            //    throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            //foreach (var m in analyzer.FactorizedMatrices.Values)
            //    m.Solve(vIn, ((Vector<double>)vOut).Data);
            preconditioner.Solve(vIn, ((Vector<double>)vOut).Data);
        }

        public void MultiplyWithMatrix(IVector<double> vIn, IVector<double> vOut)
        {
            ((SkylineMatrix2D<double>)solver.SubdomainsDictionary.Values.First().Matrix).Multiply(vIn, ((Vector<double>)vOut).Data, 1.0, 0, 0, true);
        }

        public double InitializeAndGetResidual(IList<ISolverSubdomain> subdomains, double[] r, double[] x)
        {
            double detf = 0;
            double temp = 0;

            if (subdomains.Count != 1)
                throw new InvalidOperationException("Cannot initialize and calculate residuals with more than one subdomains");

            foreach (ISolverSubdomain subdomain in subdomains)
            {
                Array.Copy(((Vector<double>)subdomain.RHS).Data, r, subdomain.RHS.Length);
                //subdomain.SubdomainToGlobalVector(((Vector<double>)subdomain.RHS).Data, r);
                var s = (SkylineMatrix2D<double>)subdomain.Matrix;

                if (preconditioner == null)
                {
                    preconditioner = (SkylineMatrix2D<double>)s.Clone();
                    preconditioner.Factorize(1e-8, new List<Vector<double>>(), new List<int>());
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
