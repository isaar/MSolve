using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Discretization
{
	public class GenericDOFEnumerator : IElementDOFEnumerator
	{
		public IList<IList<DOFType>> GetDOFTypes(IElement element)
		{
			return element.IElementType.GetElementDOFTypes(element);
		}

		public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement element)
		{
			return element.IElementType.GetElementDOFTypes(element);
		}

		public IList<INode> GetNodesForMatrixAssembly(IElement element)
		{
			return element.INodes;
		}

		public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
		{
			return matrix;
		}

		public double[] GetTransformedDisplacementsVector(double[] vector)
		{
			return vector;
		}

        public double[] GetTransformedForcesVector(double[] vector)
        { 
              return vector; 
        }

    }
}
