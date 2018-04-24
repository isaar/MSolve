using System;
using System.Collections.Generic;
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
        public SkylineSolver(Model2D model) : base(model) { }

        public override void Solve()
        {
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            var assembler = new GlobalSkylineAssembler();
            (SkylineMatrix Kuu, CSRMatrix Kuc) = assembler.BuildGlobalMatrix(model, DOFEnumerator);
            Vector rhs = CalcEffectiveRhs(Kuc);
            Solution = Kuu.FactorCholesky(true).SolveLinearSystem(rhs);
        }
    }
}
