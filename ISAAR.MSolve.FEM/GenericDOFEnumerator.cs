using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM
{
    public class GenericDOFEnumerator : IFiniteElementDOFEnumerator
    {
        public IList<IList<DOFType>> GetDOFTypes(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
        {
            return matrix;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public double[] GetTransformedVector(double[] vector)
        {
            return vector;
        }
    }
}
