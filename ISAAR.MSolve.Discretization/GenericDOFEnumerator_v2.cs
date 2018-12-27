using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Discretization
{
	public class GenericDofEnumerator_v2 : IElementDofEnumerator_v2
	{
		public IList<IList<DOFType>> GetDOFTypes(IElement_v2 element) 
            => element.ElementType.GetElementDOFTypes(element);

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement_v2 element) 
            => element.ElementType.GetElementDOFTypes(element);

        public IList<INode> GetNodesForMatrixAssembly(IElement_v2 element) => element.Nodes;

        public IMatrix GetTransformedMatrix(IMatrix matrix) => matrix;

        public double[] GetTransformedDisplacementsVector(double[] vector) => vector;

        public double[] GetTransformedForcesVector(double[] vector) => vector;
    }
}
