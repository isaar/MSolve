using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: It would be nice to have a GetOriginalMatrix() method. Investigate what changes that needs and whether a class that 
//      does not implement IElementDofEnumerator is better suited.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis
{
    public class ScalingDofEnumerator : IElementDofEnumerator
    {
        public ScalingDofEnumerator(double scaleFactor = 1.0)
        {
            this.ScaleFactor = scaleFactor;
        }

        public double ScaleFactor { get; set; }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForMatrixAssembly(IElement element)
                    => element.ElementType.GetElementDofTypes(element);

        public IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForDofEnumeration(IElement element)
            => element.ElementType.GetElementDofTypes(element);

        public IReadOnlyList<INode> GetNodesForMatrixAssembly(IElement element) => element.Nodes;

        public IMatrix GetTransformedMatrix(IMatrix matrix)
        {
            //TODO: Returning (and returning) the same matrix is misleading. The client expects to get a new object that can be 
            //      mutated independently from the matrix passed in. On the other hand it is much faster and the default 
            //      IElementDofEnumerator does not copy anything.
            Debug.Assert(ScaleFactor != 0.0, "Zero element stiffness matrix found. This may lead to a singular global matrix.");
            if (ScaleFactor != 1.0) matrix.ScaleIntoThis(ScaleFactor); // The 1.0 case is possible when calculating the unit stiffness.
            return matrix;
        }

        public double[] GetTransformedDisplacementsVector(double[] vector) => vector; //TODO: should this be scaled?

        public double[] GetTransformedForcesVector(double[] vector) => vector.Scale(ScaleFactor); //TODO: should this be scaled?

    }
}
