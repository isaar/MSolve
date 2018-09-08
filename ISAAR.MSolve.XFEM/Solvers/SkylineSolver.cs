using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class SkylineSolver: SolverBase
    {
        public SkylineSolver(Model2D model) : base(model)
        { }

        public override void Solve()
        {
            // Interleaved dof enumerator seems to be faster. I expect it to result in reduced bandwidth compared to separate
            // dof enumerator.
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);
            var assembler = new GlobalSkylineAssembler();
            (SkylineMatrix Kuu, DokRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            Solution = Kuu.FactorCholesky(true).SolveLinearSystem(rhs);
        }
    }
}
