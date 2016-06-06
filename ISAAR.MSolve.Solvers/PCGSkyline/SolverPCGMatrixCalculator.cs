using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.Solvers.PCGSkyline
{
    public class SolverPCGMatrixCalculator<T> : ISolverPCGMatrixCalculator where T : IMatrix2D<double>
    {
        #region ISolverPCGMatrixCalculator Members

        private readonly SolverPCG<T> solver;

        public SolverPCGMatrixCalculator(SolverPCG<T> solver)
        {
            this.solver = solver;
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
            solver.Precondition(vIn, vOut);
        }

        public void MultiplyWithMatrix(IVector<double> vIn, IVector<double> vOut)
        {
            solver.MultiplyWithMatrix(vIn, vOut);
        }

        #endregion
    }
}
