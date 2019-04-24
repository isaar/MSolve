using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Entities
{
    public class XSubdomain : ISubdomain
    {
        private readonly List<XNode> nodes = new List<XNode>();

        public XSubdomain(int id)
        {
            this.ID = id;
        }

        public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

        public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get ; set; }
        public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

        IReadOnlyList<IElement> ISubdomain.Elements => Elements;
        public List<XElement> Elements { get; } = new List<XElement>();

        public Vector Forces { get; set; } //TODO: this doesn't belong here

        public int ID { get; }

        public bool MaterialsModified { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        IReadOnlyList<INode> ISubdomain.Nodes => nodes;
        public IReadOnlyList<XNode> Nodes => nodes;

        public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)
        {
            var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public double[] CalculateElementDisplacements(XElement element, IVectorView globalDisplacementVector)
        {
            double[] elementNodalDisplacements = 
                FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses() => throw new NotImplementedException();

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            throw new NotImplementedException();
        }

        public void ResetMaterialsModifiedProperty() => throw new NotImplementedException();

        public void SaveMaterialState() => throw new NotImplementedException();

        public void ScaleConstraints(double scalingFactor) => throw new NotImplementedException();
    }
}
