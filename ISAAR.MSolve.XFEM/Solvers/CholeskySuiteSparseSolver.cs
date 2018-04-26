using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class CholeskySuiteSparseSolver: SolverBase
    {
        public CholeskySuiteSparseSolver() { }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            var assembler = new GlobalDOKAssembler();
            (DOKSymmetricColMajor Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);
            using (CholeskySuiteSparse factorization = Kuu.BuildSymmetricCSCMatrix(true).FactorCholesky())
            {
                Solution = factorization.SolveLinearSystem(rhs);
            }

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
