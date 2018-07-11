using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Interfaces
{
    public interface IIsogeometricDOFEnumerator
    {
        IList<IList<DOFType>> GetDOFTypes(Element element);
        IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element);
        IList<ControlPoint> GetNodesForMatrixAssembly(Element element);
        IMatrix2D GetTransformedMatrix(IMatrix2D matrix);
        double[] GetTransformedVector(double[] vector);
    }
}
