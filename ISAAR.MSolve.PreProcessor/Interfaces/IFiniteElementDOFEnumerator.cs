using System;
using System.Collections.Generic;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IFiniteElementDOFEnumerator
    {
        IList<IList<DOFType>> GetDOFTypes(Element element);
        IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element);
        IList<Node> GetNodesForMatrixAssembly(Element element);
        IMatrix2D<double> GetTransformedMatrix(IMatrix2D<double> matrix);
        double[] GetTransformedVector(double[] vector);
    }
}
