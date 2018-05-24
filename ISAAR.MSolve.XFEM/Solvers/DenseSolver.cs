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
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class DenseSolver : SolverBase
    {
        public DenseSolver(Model2D model) : base(model)
        { }

        public override void Solve()
        {
            //DofOrderer = DofOrdererSeparate.Create(model);
            DofOrderer = InterleavedDofOrderer.Create(model);
            var assembler = new GlobalSkylineAssembler();
            (Matrix Kuu, Matrix Kuc) = GlobalDenseAssembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);
            Solution = Kuu.FactorCholesky().SolveLinearSystem(rhs);
        }
    }
}
