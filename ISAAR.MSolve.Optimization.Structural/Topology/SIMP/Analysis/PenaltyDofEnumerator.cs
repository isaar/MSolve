using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis
{
    public class PenaltyDofEnumerator : IElementDofEnumerator_v2
    {
        public double Penalty { get; set; } = double.NaN;

        public IList<IList<DOFType>> GetDOFTypes(IElement_v2 element)
                    => element.ElementType.GetElementDOFTypes(element);

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement_v2 element)
            => element.ElementType.GetElementDOFTypes(element);

        public IList<INode> GetNodesForMatrixAssembly(IElement_v2 element) => element.Nodes;

        public IMatrix GetTransformedMatrix(IMatrix matrix) => matrix.Scale(Penalty); //TODO: it would be faster to not create a new object.

        public double[] GetTransformedDisplacementsVector(double[] vector) => vector; //TODO: should this be scaled?

        public double[] GetTransformedForcesVector(double[] vector) => vector.Scale(Penalty); //TODO: should this be scaled?

    }
}
