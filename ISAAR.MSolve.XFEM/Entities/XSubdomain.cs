using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;

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
        public List<IXFiniteElement> Elements { get; } = new List<IXFiniteElement>();

        public Vector Forces { get; set; } //TODO: this doesn't belong here

        public int ID { get; }

        public bool StiffnessModified { get; set; } = true; // At first it is modified
        public bool ConnectivityModified { get; set; } = true; // At first it is modified

        IReadOnlyList<INode> ISubdomain.Nodes => nodes;
        public IReadOnlyList<XNode> Nodes => nodes;

        public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)
        {
            var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public double[] CalculateElementDisplacements(IXFiniteElement element, IVectorView globalDisplacementVector)
        {
            double[] elementNodalDisplacements = 
                FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses() => throw new NotImplementedException();

        public void DefineNodesFromElements()
        {
            nodes.Clear();
            var nodeComparer = Comparer<XNode>.Create((XNode node1, XNode node2) => node1.ID - node2.ID);
            var nodeSet = new SortedSet<XNode>(nodeComparer);
            foreach (IXFiniteElement element in Elements)
            {
                foreach (XNode node in element.Nodes) nodeSet.Add(node);
            }
            nodes.AddRange(nodeSet);
        }

        public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
        {
            //TODO: perhaps it is more efficient to traverse the global constraints instead of the subdomain's nodes, provided
            //      the latter are stored as a set. 
            //TODO: the next could be a Table method: Table.KeepDataOfRows(IEnumerable<TRow> rows)
            foreach (XNode node in Nodes)
            {
                bool isNodeConstrained = globalConstraints.TryGetDataOfRow(node,
                    out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
                if (isNodeConstrained)
                {
                    foreach (var dofDisplacementPair in constraintsOfNode)
                    {
                        Constraints[node, dofDisplacementPair.Key] = dofDisplacementPair.Value;
                    }
                }
            }
        }

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            throw new NotImplementedException();
        }

        public void ResetMaterialsModifiedProperty() => throw new NotImplementedException();

        public void SaveMaterialState() => throw new NotImplementedException();

        public void ScaleConstraints(double scalingFactor) => throw new NotImplementedException();
    }
}
