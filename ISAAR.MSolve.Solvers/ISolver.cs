using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: perhaps the solver should expose the assembler, instead of wrapping it. The assembler's interface would have to be 
//      simplified a bit though. That would violate LoD, but so does MSolve in general.
namespace ISAAR.MSolve.Solvers
{
    /// <summary>
    /// Helps translate the physical model into a linear system and then solves the latter. 
    /// </summary>
    public interface ISolver : ISystemMatrixObserver
    {
        /// <summary>
        /// A dictionary that maps subdomain ids to linear systems.
        /// </summary>
        IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }

        /// <summary>
        /// Logs information, such as linear system size, the time required for various solver tasks, etc.
        /// </summary>
        SolverLogger Logger { get; }

        /// <summary>
        /// The name of the solver for logging purposes.
        /// </summary>
        string Name { get; }

        /// <summary>
        /// Assembles the matrix that corresponds to the free freedom degrees of each whole subdomain from the matrices of its 
        /// elements.
        /// </summary>
        /// <param name="elementMatrixProvider">
        /// Determines the matrix calculated for each element (e.g. stiffness, mass, etc.)
        /// </param>
        Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider);

        /// <summary>
        /// Assembles the matrices that correspond to the free and constrained freedom degrees of each whole subdomain 
        /// from the matrices of its elements. If we denote the matrix as A, the free dofs as f and the constrained dofs as c
        /// then: A = [ Aff Acf^T; Acf Acc ] (Matlab notation). This method returns Aff, Afc, Acf, Acc. If the linear system is 
        /// symmetric, then Afc = Acf^T. In this case, these entries are only stored once and shared between the returned 
        /// Afc, Acf.
        /// </summary>
        /// <param name="elementMatrixProvider">
        /// Determines the matrix calculated for each element (e.g. stiffness, mass, etc.)
        /// </param>
        Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree, 
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider);

        /// <summary>
        /// Distributes the nodal loads defined by the preprocessor to each subdomain.
        /// </summary>
        /// <param name="globalNodalLoads">
        /// A collection of loads applied to nodes along certain freedom degrees. These are usually defined by the pre-processor 
        /// or other entities of the analysis.
        /// </param>
        Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, IDofType, double> globalNodalLoads);

        /// <summary>
        /// Initializes the state of this <see cref="ISolver"/> instance. This needs to be called only once, since it  
        /// could potentially perform actions that must not be repeated or are too expensive.
        /// </summary>
        void Initialize();

        /// <summary>
        /// Solves multiple linear systems A * X = B, where: A is one of the matrices stored in <see cref="LinearSystems"/>,
        /// B is the corresponding matrix in <paramref name="otherMatrix"/> and X is the corresponding matrix that will be 
        /// calculated as the result of inv(A) * B. 
        /// </summary>
        /// <param name="otherMatrix">
        /// The right hand side matrix for each subdomain. If the linear systems are A * X = B, then B is one of the matrices in
        /// <paramref name="otherMatrix"/>.</param>
        Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix);

        /// <summary>
        /// Orders the free and optionally the constrained freedom degrees of the model. Also remember to reset the linear 
        /// systems. 
        /// </summary>
        /// <param name="alsoOrderConstrainedDofs">
        /// If true, the constrained dofs will also be ordered. Else, <see cref="ISubdomain.ConstrainedDofOrdering"/> will
        /// be null.
        /// </param>
        void OrderDofs(bool alsoOrderConstrainedDofs);
        //TODO: Would it be better if the solver didn't modify the model? It would return the ordering 
        //      and the analyzer would. However the assembler should still be notified.
        //TODO: I think this should also reset the linear systems. Is there any chance that OrderDofs() should be called but 
        //      ILinearSystem.Reset() shouldn't?

        /// <summary>
        /// Notifies this <see cref="ISolver"/> that it cannot overwrite the data of <see cref="ILinearSystem.Matrix"/>.
        /// Some solvers would otherwise overwrite the matrices (e.g. with the factorization) to avoid using extra memory.
        /// </summary>
        void PreventFromOverwrittingSystemMatrices();

        /// <summary>
        /// Solves the linear systems.
        /// </summary>
        void Solve();
    }
}
