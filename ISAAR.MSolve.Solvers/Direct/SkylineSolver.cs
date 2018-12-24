using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: factorization should be done during Solve() to time correctly. Alternatively, an observer should record the durations.
//TODO: investigate if it is possible to avoid casting the matrix provided by the analyzer/assembler into skyline. Perhaps the
//      the matrix could be obtained directly from the assembler and the analyzer could provide delegates with the operations it 
//      wants done on the matrix, instead of doing them itself.
//TODO: try to abstract the subdomain logic from ther analyzers. I doubt it is possible though.
//TODO: directly pass the single linear system instead of a list that must be checked. The same holds for all solvers and 
//      assemblers.
namespace ISAAR.MSolve.Solvers.Direct
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on sparse symmetric positive definite 
    /// matrices stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver : ISolver_v2
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly IStructuralModel_v2 model;
        private readonly ISubdomain_v2 subdomain;
        private readonly double factorizationPivotTolerance;
        private readonly SkylineSystem linearSystem;

        private bool mustFactorize = true;
        private CholeskySkyline factorizedMatrix;

        private SkylineSolver(IStructuralModel_v2 model, double factorizationPivotTolerance, IDofOrderer dofOrderer)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            subdomain = model.Subdomains[0];

            linearSystem = new SkylineSystem(subdomain);
            LinearSystems = new Dictionary<int, ILinearSystem_v2>() { { subdomain.ID, linearSystem } };
            linearSystem.MatrixObservers.Add(this);

            this.factorizationPivotTolerance = factorizationPivotTolerance;
            this.DofOrderer = dofOrderer;
        }

        public IDofOrderer DofOrderer { get; }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        public void Initialize() { }

        public void OnMatrixSetting()
        {
            mustFactorize = true;
            factorizedMatrix = null;
        }

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public void Solve()
        {
            //TODO: This should be handled by the linear system when the dof ordering changes. Here we should just call 
            //      Solution.Clear() (or nothing at all for release builds, Solution.Clear() for debug builds.).
            if (linearSystem.Solution == null) linearSystem.Solution = linearSystem.CreateZeroVector();
            else if (HaveSubdomainDofsChanged()) linearSystem.Solution = linearSystem.CreateZeroVector();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this

            if (mustFactorize)
            {
                factorizedMatrix = linearSystem.Matrix.FactorCholesky(true, factorizationPivotTolerance); 
                mustFactorize = false;
                linearSystem.IsMatrixOverwrittenBySolver = true;
            }

            factorizedMatrix.SolveLinearSystem(linearSystem.RhsVector, linearSystem.Solution);
        }

        //TODO: Create a method in Subdomain (or its DofOrderer) that exposes whether the dofs have changed.
        /// <summary>
        /// The number of dofs might have been changed since the previous Solution vector had been created.
        /// </summary>
        private bool HaveSubdomainDofsChanged() => subdomain.DofOrdering.NumFreeDofs == linearSystem.Solution.Length;

        //TODO: Copied from Stavroulakis code. Find out what the purpose of this is. I suspect he wanted to compare with some 
        //      old solution that used single precision.
        //private void DestroyAccuracy(ILinearSystem_v2 linearSystem)
        //{
        //    if (AccuracyDigits < 1) return;

        //    for (int i = 0; i < linearSystem.RhsVector.Length; i++)
        //    {
        //        //ScientificDouble s = ScientificDouble.GetScientificDouble(subdomain.RHS[i]);
        //        //s.ReduceAccuracy(AccuracyDigits);
        //        //linearSystem.RhsVector[i] = ScientificDouble.GetDouble(s);
        //        linearSystem.RhsVector[i] = Double.Parse(String.Format("{0:" + stringFormat + "}", linearSystem.RhsVector[i]));
        //    }
        //}

        public class Builder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } 
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            public SkylineSolver BuildSolver(IStructuralModel_v2 model)
            {
                return new SkylineSolver(model, FactorizationPivotTolerance, DofOrderer);
            }
        }

        private class SkylineSystem : LinearSystem_v2<SkylineMatrix, Vector>
        {
            internal SkylineSystem(ISubdomain_v2 subdomain) : base(subdomain) { }
            public override Vector CreateZeroVector() => Vector.CreateZero(Subdomain.DofOrdering.NumFreeDofs);
            public override void GetRhsFromSubdomain() => RhsVector = Subdomain.Forces;
        }
    }
}
