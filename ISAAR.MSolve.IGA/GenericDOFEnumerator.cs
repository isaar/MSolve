using ISAAR.MSolve.IGA.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.IGA
{
    public class GenericDOFEnumerator : IIsogeometricDOFEnumerator
    {
        public IList<IList<DOFType>> GetDOFTypes(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IList<ControlPoint> GetNodesForMatrixAssembly(Element element)
        {
            return element.ControlPoints;
        }

        public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
        {
            return matrix;
        }

        public double[] GetTransformedVector(double[] vector)
        {
            return vector;
        }
    }
}
