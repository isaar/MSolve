using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public abstract class SingleSubdomainRectangularSolverBase<TMatrix>:ISolver_v2
	where TMatrix: class,IMatrix
	{
		private IGlobalMatrixRectangularAssembler<TMatrix> assembler;
		private IStructuralAsymmetricModel model;
		private ISubdomain_v2 subdomain;
		private SingleSubdomainSystem<TMatrix> linearSystem;
		private IAsymmetricDofOrderer dofRowOrderer;
		private IDofOrderer dofColOrderer;

		protected SingleSubdomainRectangularSolverBase(IStructuralAsymmetricModel model, IAsymmetricDofOrderer dofRowOrderer,
			IDofOrderer dofColOrderer, IGlobalMatrixRectangularAssembler<TMatrix> assembler, string name)
		{
			if (model.Subdomains.Count != 1) throw new InvalidSolverException(
				$"{name} can be used if there is only 1 subdomain");
			this.model = model;
			subdomain= (ISubdomain_v2)model.Subdomains[0];
			linearSystem= new SingleSubdomainSystem<TMatrix>(subdomain);
			LinearSystems= new Dictionary<int, ILinearSystem_v2>(){{subdomain.ID, linearSystem}};
			linearSystem.MatrixObservers.Add(this);

			this.dofRowOrderer = dofRowOrderer;
			this.dofColOrderer = dofColOrderer;
			this.assembler = assembler;
		}

		public abstract void HandleMatrixWillBeSet();

		public IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

		public IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider)
			=> assembler.BuildGlobalMatrix(((IAsymmetricSubdomain)subdomain).DofRowOrdering, ((IAsymmetricSubdomain)subdomain).DofColOrdering, subdomain.Elements,
				elementMatrixProvider);

		public abstract void Initialize();

		public void OrderDofsAndClearLinearSystems()
		{
			IGlobalFreeDofOrdering globalRowOrdering = dofRowOrderer.OrderDofs(model);
			IGlobalFreeDofOrdering globalColOrdering = dofColOrderer.OrderDofs(model);
			assembler.HandleDofOrderingWillBeModified();
			HandleMatrixWillBeSet();
			linearSystem.Clear();
			linearSystem.Size = globalColOrdering.SubdomainDofOrderings[subdomain].NumFreeDofs;

			model.GlobalRowDofOrdering = globalRowOrdering;
			model.GlobalColDofOrdering = globalColOrdering;

			foreach (var subdomain in model.Subdomains)
			{
				subdomain.DofRowOrdering = globalRowOrdering.SubdomainDofOrderings[(ISubdomain_v2)subdomain];
				subdomain.DofColOrdering = globalColOrdering.SubdomainDofOrderings[(ISubdomain_v2)subdomain];
			}
		}

		public void ResetSubdomainForcesVector()
		{
			foreach (ISubdomain_v2 subdomain in model.Subdomains)
			{
				subdomain.Forces = linearSystem.CreateZeroVector(); 
			}
		}
		public abstract void PreventFromOverwrittingSystemMatrices();
		public abstract void Solve();
	}
}
