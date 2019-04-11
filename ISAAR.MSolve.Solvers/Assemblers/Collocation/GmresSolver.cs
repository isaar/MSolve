using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
    public class GmresSolver : SingleSubdomainRectangularSolverBase<CsrMatrix>
    {

        public GmresSolver(IStructuralAsymmetricModel model,  AsymmetricDofOrderer dofRowOrderer, IDofOrderer dofColOrderer,
            IGlobalMatrixRectangularAssembler<CsrMatrix> assembler, string name) : base(model, dofRowOrderer, dofColOrderer, new CsrRectangularAssembler(true), "GmresSolver")
        {
        }

        public override void Initialize() { }

        public override void HandleMatrixWillBeSet()
        {
        }

        public override void PreventFromOverwrittingSystemMatrices()
        {
        }

        public override void Solve()
        {
            throw new NotImplementedException();
        }
    }
}