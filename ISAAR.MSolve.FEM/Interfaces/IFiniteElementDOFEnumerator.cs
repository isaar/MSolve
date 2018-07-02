using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElementDOFEnumerator
    {
        IList<IList<DOFType>> GetDOFTypes(Element element);
        IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element);
        IList<Node> GetNodesForMatrixAssembly(Element element);
        IMatrix2D GetTransformedMatrix(IMatrix2D matrix);
        double[] GetTransformedVector(double[] vector);
    }
}
