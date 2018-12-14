using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
	public interface IElementDOFEnumerator
	{
        /// <summary>
        /// These are the dofs of the nodes returned by <see cref="GetNodesForMatrixAssembly"/>
        /// </summary>
		IList<IList<DOFType>> GetDOFTypes(IElement element);

        /// <summary>
        /// The returned outer list will include nested lists for all <see cref="IElement.INodes"/>. When using embedding, the
        /// nested lists, that correspond to embedded nodes, will be empty.
        /// </summary>
		IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement element);

        /// <summary>
        /// When using embedding, these are the nodes of the superelement: nodes that have not been embedded and (right now all) 
        /// nodes of the host element. 
        /// </summary>
		IList<INode> GetNodesForMatrixAssembly(IElement element);

        /// <summary>
        /// Kembedded = transpose(T) * Koriginal * T
        /// </summary>
		IMatrix2D GetTransformedMatrix(IMatrix2D matrix);

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