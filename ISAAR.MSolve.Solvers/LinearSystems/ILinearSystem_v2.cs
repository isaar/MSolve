using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: setters should be removed, since they force the implementation (or the solver) to cast. Currently there is no
//      way around it though.
//TODO: This should hide distributed systems arising in domain decomposition methods. That should be transparent to the analyzers
//      and providers.
//TODO: Perhaps the matrix should not necessarily be IMatrix. It could also be a DOK. However that might complicate the
//      analyzers, providers.
#region design comments
// How this works: 
// 1) The assemblers and solvers have full access to the concrete matrices and vectors.
// The analyzers and providers can freely call non-mutating methods on the matrix and vector interface. They can also call
// mutating methods. 
// 2) If a class mutates the matrix or vectors, then the LinearSystem object must handle its state and 
// notify the other classes (solvers, assemblers, analyzers, providers) that compose it. It may also need to ask permission before
// calling a mutating method (e.g. the solver might mark the matrix for clearing, but the analyzer might still need it, 
// so the clearing request will not go through). 
// 3) Instead of using flags to determine the state of the matrix and vectors, it would be better to have solver, analyzer 
// and provider observe the LinearSystem. When someone alters e.g. the matrix, all the others will be notified to manage 
// their own memory. This system also allows communicating and approving the clearing requests. However it would be very 
// inefficient to alert all observers everytime a setter is called, therefore access to mutating methods should be controlled.
// An alternative would be to use the State pattern.
// 4) To keep analyzers and providers non-generic, LinearSystem does not provide setters that operate on the interface vectors. 
// Instead it must provide methods to initialize zero vectors (or scalar * ones()). Currently the vectors created inside the 
// analyzer and provider are all of type Vector. However that is restricting, inefficient and will cause problems once 
// distributed computing is introduced.
// 5) The system matrix is built by the assembler, thus analyzer/providers delegate setting it to the assembler. 
// 6) After initialing the matrices and vectors (setting them is done by the LinearSystem or assembler), the analyzers may call
// mutating methods on them, as described above.
// 7) What about auxiliary global matrices in analyzers/providers: these will have a trivial pattern (diagonal, dense) or  
// the same one as the system matrix. Creating them should also be delegated to the assembler (at least in the second case). 
// I am not sure yet if the analyzers/providers should have any control on the pattern (e.g. by telling the assembler to use 
// the same pattern, without specifying which one, for mass matrix as for stiffness matrix).
// 8) What about auxiliary global vectors in analyzers/providers? Instantiation should be done with abstract factory pattern 
// (LinearSystem could be the factory itself, but that violates SRP). The vector length could also be accessed from there. 
// Distributing them to and assembling them from subdomains should be hidden from the analyzers/providers as well.
// 9) An alternative to prohibition of setting the rhs vector would be to check the validity of the new value in the setters.
// As an example, the linear system initializes a vector f, giving it a GUID. The analyzer performs various operations on it
// and then sets linearSystem.RhsVector = f. Then the setter will check the GUID, cast the vector, set the concrete object,
// and notify observers that a value has been set. Instead of a GUID, the linear system can keep track of the vectors it has
// created and use that registry to accept setting its rhs vector. A similar logic can be used for operations between
// matrices with the same sparsity pattern.
#endregion
namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// A system of linear equations. It consists of a square matrix, a right hand side vector and a solution (or left 
    /// hand side) vector. In general objects implementing IAnalyzer determine the matrix and right hand side vector,  
    /// while objects implementing ISolver calculate the solution vector.
    /// </summary>
    public interface ILinearSystem_v2
    {
        /// <summary>
        /// The linear system matrix.
        /// </summary>
        IMatrixView Matrix { get; set; }

        /// <summary>
        /// Objects that will be notified when <see cref="Matrix"/> is modified.
        /// </summary>
        HashSet<ISystemMatrixObserver> MatrixObservers { get; }

        /// <summary>
        /// The right hand side vector.
        /// </summary>
        IVector RhsVector { get; set; }

        /// <summary>
        /// The number of linear equations, which is equal to the number of unknowns and the dimensions of the matrix 
        /// and vectors.
        /// </summary>
        int Size { get; }

        /// <summary>
        /// The subdomain corresponding to this linear system.
        /// </summary>
        ISubdomain_v2 Subdomain { get; } //TODO: delete this once subdomains have been abstracted.

        /// <summary>
        /// The solution vector.
        /// </summary>
        IVectorView Solution { get; }

        /// <summary>
        /// Initializes a new vector with zero entries. Its pattern depends on the solver used. The freedom degrees must be 
        /// ordered before this method can be called. Special attention is needed if the freedom degrees change during 
        /// the analysis (e.g. adaptive FEM, XFEM). 
        /// </summary>
        IVector CreateZeroVector();

        /// <summary>
        /// Clears all data stored in this <see cref="ILinearSystem_v2"/> instance and sets its <see cref="Size"/> equal to the
        /// number of freedom degrees of the corresponding <see cref="Subdomain"/>.
        /// </summary>
        void Reset();
    }
}
