using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    public class MonteCarloSolverPCGDirectMatrixCalculator_v2
    {
        private MonteCarloAnalyzerWithStochasticMaterial analyzer;
        private LdlSkyline preconditioner;

        public int VectorSize => LinearSystem.RhsVector.Length;

        public ILinearSystem_v2 LinearSystem { get; set; }

        public void Precondition(Vector rhs, Vector lhs)
        {
            //if (analyzer.FactorizedMatrices.Count != 1)
            //    throw new InvalidOperationException("Cannot precondition with more than one subdomains");

            //foreach (var m in analyzer.FactorizedMatrices.Values)
            //    m.Solve(vIn, ((Vector<double>)vOut).Data);
            preconditioner.SolveLinearSystem(rhs, lhs);
        }

        public void MultiplyWithMatrix(IVectorView lhs, IVector rhs)
        {
            //((SkylineMatrix2D)solver.SubdomainsDictionary.Values.First().Matrix).Multiply(vIn, ((Vector)vOut).Data, 1.0, 0, 0, true);
            LinearSystem.Matrix.MultiplyIntoResult(lhs, rhs);
        }

        public double InitializeAndGetResidual(IList<ILinearSystem_v2> linearSYstems, IVector res, IVector x0)
        {
            if (linearSYstems.Count != 1)
                throw new InvalidOperationException("Cannot initialize and calculate residuals with more than one subdomains");

            foreach (ILinearSystem_v2 linearSystem in linearSYstems)
            {
                res.CopyFrom(linearSystem.RhsVector);
                //Array.Copy(((Vector)subdomain.RHS).Data, r, subdomain.RHS.Length);

                //subdomain.SubdomainToGlobalVector(((Vector<double>)subdomain.RHS).Data, r);

                if (preconditioner == null)
                {
                    preconditioner = ((SkylineMatrix)linearSystem.Matrix).FactorLdl(false, 1e-8);
                }
            }

            return res.Norm2();
        }
    }
}
