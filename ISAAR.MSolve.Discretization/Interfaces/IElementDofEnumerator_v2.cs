using ISAAR.MSolve.LinearAlgebra.Matrices;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IElementDofEnumerator_v2
	{
        /// <summary>
        /// These are the dofs of the nodes returned by <see cref="GetNodesForMatrixAssembly"/>
        /// </summary>
		IList<IList<DOFType>> GetDOFTypes(IElement_v2 element); //TODO: rename GetDofTypesForMatrixAssembly()

        /// <summary>
        /// The returned outer list will include nested lists for all <see cref="IElement_v2.Nodes"/>. When using embedding, the
        /// nested lists, that correspond to embedded nodes, will be empty.
        /// </summary>
		IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement_v2 element);

        /// <summary>
        /// When using embedding, these are the nodes of the superelement: nodes that have not been embedded and (right now all) 
        /// nodes of the host element. 
        /// </summary>
		IList<INode> GetNodesForMatrixAssembly(IElement_v2 element);

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