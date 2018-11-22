using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using IVectorOLD = ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces.IVector;
using VectorOLD = ISAAR.MSolve.Numerical.LinearAlgebra.Vector;

//TODO: remove code that calculates rhs vector components. It should be moved to dedicated classes like EquivalentLoadAssembler,
//      so that it can be reused between subdomains of different projects (FEM, IGA, XFEM).
namespace ISAAR.MSolve.FEM.Entities
{
    public class Subdomain_v2 : ISubdomain_v2
    {
        private readonly List<Node> nodes = new List<Node>();

        public Subdomain_v2(int id)
        {
            this.ID = id;
        }

        public Table<INode, DOFType, double> Constraints { get; } = new Table<INode, DOFType, double>();

        IReadOnlyList<IElement> ISubdomain_v2.Elements => Elements;
        public List<Element> Elements { get; } = new List<Element>();

        //public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

        public int ID { get; }

        IReadOnlyList<INode> ISubdomain_v2.Nodes => nodes;
        public IReadOnlyList<Node> Nodes => nodes;

        public ISubdomainFreeDofOrdering DofOrdering { get; set; }

        public double[] Forces { get; set; } //TODO: this doesn't belong here

        public bool MaterialsModified { get; set; }
        //public bool MaterialsModified
        //{
        //    get
        //    {
        //        bool modified = false;
        //        foreach (Element element in elementsDictionary.Values)
        //            if (element.ElementType.MaterialModified)
        //            {
        //                modified = true;
        //                break;
        //            }
        //        return modified;
        //    }
        //}

        //TODO: Ideally this is set by the Model, Cluster and should not be modified during the analysis.
        public Table<Node, DOFType, double> NodalLoads { get; set; }

        public void BuildConstraintDisplacementDictionary()
        {
            //TODO: constraints should not be saved inside the nodes, just have a Table stored by Model.
            foreach (Node node in Nodes)
            {
                if (node.Constraints == null) continue;
                foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            }
        }

        public void BuildNodesList()
        {
            var nodeComparer = Comparer<Node>.Create((Node node1, Node node2) => node1.ID - node2.ID);
            var nodeSet = new SortedSet<Node>(nodeComparer);
            foreach (Element element in Elements)
            {
                foreach (Node node in element.Nodes) nodeSet.Add(node);
            }
            nodes.AddRange(nodeSet);

            //List<int> nodeIDs = new List<int>();
            //Dictionary<int, Node> nodes = new Dictionary<int, Node>();
            //foreach (Element element in Elements)
            //    foreach (Node node in element.Nodes)
            //    {
            //        nodeIDs.Add(node.ID);
            //        if (!nodes.ContainsKey(node.ID))
            //            nodes.Add(node.ID, node);
            //    }

            //nodeIDs = new List<int>(nodeIDs.Distinct<int>());
            //nodeIDs.Sort();
            //foreach (int nodeID in nodeIDs)
            //    nodesDictionary.Add(nodeID, nodes[nodeID]);

            ////foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
            ////    EmbeddedNodes.Add(e);
        }

        //TODO: This belongs in EquivalentLoadsAssembler
        //TODO: the constraintScalingFactor parameter is not used.
        public Vector CalculateElementIncrementalConstraintDisplacements(Element element, double constraintScalingFactor)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            int numElementDofs = DofOrdering.CountElementDofs(element);
            var elementNodalDisplacements = Vector.CreateZero(numElementDofs);
            ApplyConstraintDisplacements(element, elementNodalDisplacements);
            return elementNodalDisplacements;
        }

        public Vector CalculateElementNodalDisplacements(Element element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            var elementNodalDisplacements = Vector.CreateZero(DofOrdering.CountElementDofs(element));
            DofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector, elementNodalDisplacements);
            ApplyConstraintDisplacements(element, elementNodalDisplacements);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses()
        {
            foreach (Element element in Elements) element.ElementType.ClearMaterialStresses();
        }

        public void EnumerateDOFs()
        {
            //DofOrdering = dofOrderer(this);
            //Forces = new double[DofOrdering.NumFreeDofs];
        }

        public int[] GetCornerNodes()
        {
            int nodex1y1z1 = -1;
            int nodex2y1z1 = -1;
            int nodex1y2z1 = -1;
            int nodex2y2z1 = -1;
            int nodex1y1z2 = -1;
            int nodex2y1z2 = -1;
            int nodex1y2z2 = -1;
            int nodex2y2z2 = -1;
            double x1 = Double.MaxValue;
            double x2 = Double.MinValue;
            double y1 = Double.MaxValue;
            double y2 = Double.MinValue;
            double z1 = Double.MaxValue;
            double z2 = Double.MinValue;

            foreach (var kv in nodes)
            {
                if (x1 > kv.X) x1 = kv.X;
                if (x2 < kv.X) x2 = kv.X;
                if (y1 > kv.Y) y1 = kv.Y;
                if (y2 < kv.Y) y2 = kv.Y;
                if (z1 > kv.Z) z1 = kv.Z;
                if (z2 < kv.Z) z2 = kv.Z;

                if (x1 == kv.X && y1 == kv.Y && z1 == kv.Z) nodex1y1z1 = kv.ID;
                if (x2 == kv.X && y1 == kv.Y && z1 == kv.Z) nodex2y1z1 = kv.ID;
                if (x1 == kv.X && y2 == kv.Y && z1 == kv.Z) nodex1y2z1 = kv.ID;
                if (x2 == kv.X && y2 == kv.Y && z1 == kv.Z) nodex2y2z1 = kv.ID;
                if (x1 == kv.X && y1 == kv.Y && z2 == kv.Z) nodex1y1z2 = kv.ID;
                if (x2 == kv.X && y1 == kv.Y && z2 == kv.Z) nodex2y1z2 = kv.ID;
                if (x1 == kv.X && y2 == kv.Y && z2 == kv.Z) nodex1y2z2 = kv.ID;
                if (x2 == kv.X && y2 == kv.Y && z2 == kv.Z) nodex2y2z2 = kv.ID;
            }

            return new[] { nodex1y1z1, nodex2y1z1, nodex1y2z1, nodex2y2z1, nodex1y1z2, nodex2y1z2, nodex1y2z2, nodex2y2z2 };
        }

        public IVector GetRHSFromSolution(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(DofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in Elements)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

                //TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
                double[] localSolution = CalculateElementNodalDisplacements(element, solution).ToRawArray();
                double[] localdSolution = CalculateElementNodalDisplacements(element, dSolution).ToRawArray();
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain_v2.MaterialsModified = true;
                var f = Vector.CreateFromArray(element.ElementType.CalculateForces(element, localSolution, localdSolution));
                DofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.MaterialsModified = false;
            foreach (Element element in Elements) element.ElementType.ResetMaterialModified();
        }

        public void SaveMaterialState()
        {
            foreach (Element element in Elements) element.ElementType.SaveMaterialState();
        }

        public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

        private void ApplyConstraintDisplacements(Element element, Vector elementNodalDisplacements)
        {
            int elementDofIdx = 0;
            IList<INode> nodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> dofs = element.ElementType.DOFEnumerator.GetDOFTypes(element);
            for (int i = 0; i < nodes.Count; ++i)
            {
                //bool isConstrainedNode = constraintsDictionary.TryGetValue(nodes[i].ID, 
                //    out Dictionary<DOFType, double> constrainedDOFs);
                bool isConstrainedNode = Constraints.TryGetDataOfRow(nodes[i],
                    out IReadOnlyDictionary<DOFType, double> constrainedDOFs);
                if (isConstrainedNode)
                {
                    foreach (DOFType dofType in dofs[i])
                    {
                        bool isConstrainedDof = constrainedDOFs.TryGetValue(dofType, out double constraintDisplacement);
                        //if (isConstrainedNode && isConstrainedDof)
                        if (isConstrainedDof)
                        {
                            Debug.Assert(elementNodalDisplacements[elementDofIdx] == 0); // TODO: and why is this an assumption?
                            elementNodalDisplacements[elementDofIdx] = constraintDisplacement;
                        }
                        ++elementDofIdx;
                    }
                }
                else elementDofIdx += dofs[i].Count;
            }
        }
    }
}
