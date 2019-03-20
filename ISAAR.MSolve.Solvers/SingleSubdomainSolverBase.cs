using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
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

        public virtual Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider_v2 elementMatrixProvider)
        {
            return new Dictionary<int, IMatrix>
            {
                { subdomain.ID,
                    assembler.BuildGlobalMatrix(subdomain.FreeDofOrdering, subdomain.Elements, elementMatrixProvider) }
            };
        }

        public virtual Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider_v2 elementMatrixProvider)
        {
            if (subdomain.ConstrainedDofOrdering == null)
            {
                throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs,"
                    + " they must have been ordered first.");
            }
            return new Dictionary<int, (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc)>
            {
                { subdomain.ID, assembler.BuildGlobalSubmatrices(subdomain.FreeDofOrdering, subdomain.ConstrainedDofOrdering,
                        subdomain.Elements, elementMatrixProvider) }
            };
            
        }

        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, DOFType, double> globalNodalLoads)
        {
            var subdomainLoads = new SortedDictionary<int, double>();
            foreach ((INode node, DOFType dofType, double amount) in globalNodalLoads)
            {
                int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                subdomainLoads[subdomainDofIdx] = amount;
            }
            return new Dictionary<int, SparseVector>
            {
                { subdomain.ID, SparseVector.CreateFromDictionary(subdomain.FreeDofOrdering.NumFreeDofs, subdomainLoads) }
            };
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            if (otherMatrix.Count != 1) throw new InvalidSolverException("There can only be 1 subdomain when using this solver");
            KeyValuePair<int, IMatrixView> idMatrixPair = otherMatrix.First();
            int id = idMatrixPair.Key;
            Debug.Assert(id == subdomain.ID, 
                "The matrix that will be multiplied with the inverse system matrix belongs to a different subdomain.");
            Matrix result = InverseSystemMatrixTimesOtherMatrix(idMatrixPair.Value);
            return new Dictionary<int, Matrix>() { { id, result } };
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
        {
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            assembler.HandleDofOrderingWillBeModified();
            
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();
        }

        public abstract void Initialize();
        public abstract void HandleMatrixWillBeSet();
        public abstract void PreventFromOverwrittingSystemMatrices();
        public abstract void Solve();
        protected abstract Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix);
    }
}
