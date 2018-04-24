using System;
using System.Collections.Generic;
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
    class DenseSolver : ISolver
    {
        private readonly Model2D model;

        public DenseSolver(Model2D model)
        {
            this.model = model;
        }

        public IDOFEnumerator DOFEnumerator { get; private set; }

        public Vector Solution { get; private set; }

        public void Initialize() { }

        public void Solve()
        {
            //DOFEnumerator = DOFEnumeratorSeparate.Create(model);
            DOFEnumerator = DOFEnumeratorInterleaved.Create(model);
            (Matrix matrix, Vector rhs) = ReduceToSimpleLinearSystem();
            Solution = matrix.FactorCholesky().SolveLinearSystem(rhs);
        }

        /// <summary>
        /// The extended linear system is:
        /// [Kcc Kcu; Kuc Kuu] * [uc; uu] = [Fc; Fu]
        /// where c are the standard constrained dofs, f are the standard free dofs, e are the enriched dofs and 
        /// u = Union(f,c) are both the dofs with unknown left hand side vectors: uu = [uf; ue].
        /// To solve the system (for the unknowns ul):
        /// i) Kuu * uu = Fu - Kuc * uc = Feff
        /// ii) uu = Kuu \ Feff 
        /// </summary>
        /// <returns></returns>
        private (Matrix matrix, Vector rhs) ReduceToSimpleLinearSystem()
        {
            (Matrix Kuu, Matrix Kuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, DOFEnumerator);

            // TODO: Perhaps a dedicated class should be responsible for these vectors
            Vector Fu = model.CalculateFreeForces(DOFEnumerator);
            Vector uc = model.CalculateConstrainedDisplacements(DOFEnumerator);
            Vector Feff = Fu - Kuc * uc;
            return (Kuu, Feff);
        }
    }
}
