using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class SkylineSolver: SolverBase
    {
        public SkylineSolver() { }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Interleaved dof enumerator seems to be faster. I expect it to result in reduced bandwidth compared to separate
            // dof enumerator.
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            var assembler = new GlobalSkylineAssembler();
            (SkylineMatrix Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);
            Solution = Kuu.FactorCholesky(true).SolveLinearSystem(rhs);

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
