using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class DenseSolver : SolverBase
    {
        public DenseSolver(){ }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            var assembler = new GlobalSkylineAssembler();
            (Matrix Kuu, Matrix Kuc) = GlobalDenseAssembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);
            Solution = Kuu.FactorCholesky().SolveLinearSystem(rhs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
