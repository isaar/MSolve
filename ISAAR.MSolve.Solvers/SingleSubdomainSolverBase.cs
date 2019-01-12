using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: not sure this class should be in this namespace
namespace ISAAR.MSolve.Solvers
{
    /// <summary>
    /// Base implementation for solver that do not use domain decomposition.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix">The type of the linear system's matrix.</typeparam>
    public abstract class SingleSubdomainSolverBase<TMatrix> : ISolver_v2
        where TMatrix : class, IMatrix
    {
        protected readonly IGlobalMatrixAssembler<TMatrix> assembler;
        protected readonly IDofOrderer dofOrderer;
        protected readonly IStructuralModel_v2 model;
        protected readonly string name; // for error messages
        protected readonly ISubdomain_v2 subdomain;
        protected readonly SingleSubdomainSystem<TMatrix> linearSystem;

        protected SingleSubdomainSolverBase(IStructuralModel_v2 model, IDofOrderer dofOrderer, 
            IGlobalMatrixAssembler<TMatrix> assembler, string name)
        {
            if (model.Subdomains.Count != 1) throw new InvalidSolverException(
                $"{name} can be used if there is only 1 subdomain");
            this.model = model;
            subdomain = model.Subdomains[0];

            linearSystem = new SingleSubdomainSystem<TMatrix>(subdomain);
            LinearSystems = new Dictionary<int, ILinearSystem_v2>() { { subdomain.ID, linearSystem } };
            linearSystem.MatrixObservers.Add(this);

            this.dofOrderer = dofOrderer;
            this.assembler = assembler;
        }

        public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        public virtual IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
            => assembler.BuildGlobalMatrix(subdomain.DofOrdering, subdomain.Elements, elementMatrixProvider);

        public virtual void OrderDofsAndClearLinearSystems()
        {
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderDofs(model);
            assembler.HandleDofOrderingWillBeModified();
            HandleMatrixWillBeSet();
            linearSystem.Clear();
            linearSystem.Size = globalOrdering.SubdomainDofOrderings[subdomain].NumFreeDofs;

            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                subdomain.DofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];

                // If we decide subdomain.Forces will always be a Vector or double[] then this process could be done elsewhere.
                subdomain.Forces = linearSystem.CreateZeroVector();
            }
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();
        }

        public abstract void Initialize();
        public abstract void HandleMatrixWillBeSet();
        public abstract void PreventFromOverwrittingSystemMatrices();
        public abstract void Solve();
    }
}
