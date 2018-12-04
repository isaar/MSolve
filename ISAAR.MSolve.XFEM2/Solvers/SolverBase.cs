using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    abstract class SolverBase: ISolver
    {
        protected readonly Model2D model;

        protected SolverBase(Model2D model)
        {
            this.model = model;
            Logger = new SolverLogger();
        }

        public IDofOrderer DofOrderer { get; protected set; }

        public SolverLogger Logger { get; }

        public IVector Solution { get; protected set; }

        public virtual void Initialize() // Many solvers do not need to initialize state
        {
        } 

        public abstract void Solve();

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
        protected Vector CalcEffectiveRhs(DokRowMajor globalUnconstrainedConstrained)
        {
            Vector Fu = model.CalculateFreeForces(DofOrderer);
            Vector uc = model.CalculateConstrainedDisplacements(DofOrderer);
            Vector Feff = Fu - globalUnconstrainedConstrained.MultiplyRight(uc);
            return Feff;
        }

        protected Vector CalcEffectiveRhs(IMatrixView globalUnconstrainedConstrained)
        {
            Vector Fu = model.CalculateFreeForces(DofOrderer);
            Vector uc = model.CalculateConstrainedDisplacements(DofOrderer);
            Fu.SubtractIntoThis(globalUnconstrainedConstrained.Multiply(uc));
            return Fu;
        }
    }
}
