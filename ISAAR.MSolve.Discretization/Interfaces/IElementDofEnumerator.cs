using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IElementDofEnumerator
	{
        /// <summary>
        /// These are the dofs of the nodes returned by <see cref="GetNodesForMatrixAssembly"/>
        /// </summary>
		IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForMatrixAssembly(IElement element);

        /// <summary>
        /// The returned outer list will include nested lists for all <see cref="IElement.Nodes"/>. When using embedding, the
        /// nested lists that correspond to embedded nodes, will be empty.
        /// </summary>
		IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForDofEnumeration(IElement element);

        /// <summary>
        /// When using embedding, these are the nodes of the superelement: nodes that have not been embedded and (right now all) 
        /// nodes of the host element. 
        /// </summary>
		IReadOnlyList<INode> GetNodesForMatrixAssembly(IElement element);

		IMatrix GetTransformedMatrix(IMatrix matrix);

        /// <summary>
        /// Returns element local displacements.
        /// </summary>
		double[] GetTransformedDisplacementsVector(double[] superElementDisplacements);

        /// <summary>
        /// Returns super-element forces.
        /// </summary>
        double[] GetTransformedForcesVector(double[] elementLocalForces);
    }
}