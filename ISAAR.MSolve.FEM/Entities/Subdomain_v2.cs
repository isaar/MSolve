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

namespace ISAAR.MSolve.FEM.Entities
{
    public class Subdomain_v2 : ISubdomain_v2
    {
        //TODO: remove these and let the solver's dof orderer do the job.
        public delegate IDofOrdering OrderDofs(ISubdomain_v2 subdomain);
        private readonly OrderDofs dofOrderer;

        //private readonly IList<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        //private readonly Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> globalNodalDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly Dictionary<int, Dictionary<DOFType, double>> constraintsDictionary = new Dictionary<int, Dictionary<DOFType, double>>();
        private double[] forces;

        public Subdomain_v2(int id, OrderDofs dofOrderer)
        {
            this.ID = id;
            this.dofOrderer = dofOrderer;
        }

        #region Properties
        //public IList<EmbeddedNode> EmbeddedNodes
        //{
        //    get { return embeddedNodes; }
        //}

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, IElement> ΙElementsDictionary
        {
            get
            {
                var a = new Dictionary<int, IElement>();
                foreach (var element in elementsDictionary.Values)
                    a.Add(element.ID, element);
                return a;
            }
        }

        //TODO: Ideally this is set by the Model, Cluster and should not be modified during the analysis.
        public Table<Node, DOFType, double> NodalLoads { get; set; }

        public Dictionary<int, Node> NodesDictionary
        {
            get { return nodesDictionary; }
        }


        public IList<Node> Nodes
        {
            get { return nodesDictionary.Values.ToList<Node>(); }
        }

        IReadOnlyList<INode> ISubdomain_v2.Nodes => nodesDictionary.Values.ToList<INode>();

        public Dictionary<int, Dictionary<DOFType, double>> Constraints => constraintsDictionary;

        public IDofOrdering DofOrdering { get; set; }

        //public Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary
        //{
        //    get { return nodalDOFsDictionary; }
        //}

        public Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary
        {
            get { return globalNodalDOFsDictionary; }
        }

        public double[] Forces
        {
            get { return forces; }
        }

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

        public int ID { get; }
        #endregion

        #region Data inteconnection routines
        public void EnumerateDOFs()
        {
            DofOrdering = dofOrderer(this);
            forces = new double[DofOrdering.NumFreeDofs];
        }

        //public void AssignGlobalNodalDOFsFromModel(Dictionary<int, Dictionary<DOFType, int>> glodalDOFsDictionary)
        //{
        //    foreach (int nodeID in nodalDOFsDictionary.Keys)
        //    {
        //        Dictionary<DOFType, int> dofTypes = nodalDOFsDictionary[nodeID];
        //        Dictionary<DOFType, int> globalDOFTypes = new Dictionary<DOFType, int>(dofTypes.Count);
        //        foreach (DOFType dofType in dofTypes.Keys)
        //            globalDOFTypes.Add(dofType, glodalDOFsDictionary[nodeID][dofType]);
        //        globalNodalDOFsDictionary.Add(nodeID, globalDOFTypes);
        //    }
        //}

        public void AssignGlobalNodalDOFsFromModel_v2(Dictionary<int, Dictionary<DOFType, int>> glodalDOFsDictionary)
        {
            foreach (Node node in DofOrdering.FreeDofs.GetRows())
            {
                IEnumerable<DOFType> subdomainDofsOfNode = DofOrdering.FreeDofs.GetColumnsOfRow(node);
                var globalDofsOfNode = new Dictionary<DOFType, int>();
                foreach (DOFType dofType in subdomainDofsOfNode)
                {
                    globalDofsOfNode.Add(dofType, glodalDOFsDictionary[node.ID][dofType]);
                }
                globalNodalDOFsDictionary.Add(node.ID, globalDofsOfNode);
            }
        }

        public void BuildNodesDictionary()
        {
            List<int> nodeIDs = new List<int>();
            Dictionary<int, Node> nodes = new Dictionary<int, Node>();
            foreach (Element element in elementsDictionary.Values)
                foreach (Node node in element.Nodes)
                {
                    nodeIDs.Add(node.ID);
                    if (!nodes.ContainsKey(node.ID))
                        nodes.Add(node.ID, node);
                }

            nodeIDs = new List<int>(nodeIDs.Distinct<int>());
            nodeIDs.Sort();
            foreach (int nodeID in nodeIDs)
                nodesDictionary.Add(nodeID, nodes[nodeID]);

            //foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
            //    embeddedNodes.Add(e);
        }

        public void BuildConstraintDisplacementDictionary()
        {
            foreach (Node node in nodesDictionary.Values)
            {
                if (node.Constraints == null) continue;
                constraintsDictionary[node.ID] = new Dictionary<DOFType, double>();
                foreach (Constraint constraint in node.Constraints)
                {
                    constraintsDictionary[node.ID][constraint.DOF] = constraint.Amount;
                }
            }
        }

        #endregion

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

            foreach (var kv in nodesDictionary)
            {
                if (x1 > kv.Value.X) x1 = kv.Value.X;
                if (x2 < kv.Value.X) x2 = kv.Value.X;
                if (y1 > kv.Value.Y) y1 = kv.Value.Y;
                if (y2 < kv.Value.Y) y2 = kv.Value.Y;
                if (z1 > kv.Value.Z) z1 = kv.Value.Z;
                if (z2 < kv.Value.Z) z2 = kv.Value.Z;

                if (x1 == kv.Value.X && y1 == kv.Value.Y && z1 == kv.Value.Z) nodex1y1z1 = kv.Key;
                if (x2 == kv.Value.X && y1 == kv.Value.Y && z1 == kv.Value.Z) nodex2y1z1 = kv.Key;
                if (x1 == kv.Value.X && y2 == kv.Value.Y && z1 == kv.Value.Z) nodex1y2z1 = kv.Key;
                if (x2 == kv.Value.X && y2 == kv.Value.Y && z1 == kv.Value.Z) nodex2y2z1 = kv.Key;
                if (x1 == kv.Value.X && y1 == kv.Value.Y && z2 == kv.Value.Z) nodex1y1z2 = kv.Key;
                if (x2 == kv.Value.X && y1 == kv.Value.Y && z2 == kv.Value.Z) nodex2y1z2 = kv.Key;
                if (x1 == kv.Value.X && y2 == kv.Value.Y && z2 == kv.Value.Z) nodex1y2z2 = kv.Key;
                if (x2 == kv.Value.X && y2 == kv.Value.Y && z2 == kv.Value.Z) nodex2y2z2 = kv.Key;
            }

            return new[] { nodex1y1z1, nodex2y1z1, nodex1y2z1, nodex2y2z1, nodex1y1z2, nodex2y1z2, nodex1y2z2, nodex2y2z2 };
        }

        public void ScaleConstraints(double scalingFactor)
        {
            var nodeIds = constraintsDictionary.Keys.ToList();
            foreach (var nodeId in nodeIds)
            {
                var dofs = constraintsDictionary[nodeId].Keys.ToList();
                foreach (DOFType dof in dofs)
                {
                    constraintsDictionary[nodeId][dof] = constraintsDictionary[nodeId][dof] * scalingFactor;
                }
            }
        }

        public double[] CalculateElementIncrementalConstraintDisplacements(Element element, double constraintScalingFactor)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            int localDOFs = 0;
            foreach (IList<DOFType> dofs in element.ElementType.DOFEnumerator.GetDOFTypes(element)) localDOFs += dofs.Count;
            var elementNodalDisplacements = new double[localDOFs];
            ApplyConstraintDisplacements(element, elementNodalDisplacements);
            var incrementalNodalDisplacements = new double[localDOFs];
            elementNodalDisplacements.CopyTo(incrementalNodalDisplacements, 0);
            
            return incrementalNodalDisplacements;
        }

        public double[] CalculateElementNodalDisplacements(Element element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            double[] elementNodalDisplacements = DofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            ApplyConstraintDisplacements(element, elementNodalDisplacements);
            return elementNodalDisplacements;
        }

        private void ApplyConstraintDisplacements(Element element, double[] elementNodalDisplacements)
        {
            int pos = 0;
            for (int i = 0; i < element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; i++)
            {
                INode node = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element)[i]; //Node node = element.Nodes[i];
                bool isConstrainedNode = constraintsDictionary.TryGetValue(node.ID, 
                    out Dictionary<DOFType, double> constrainedDOFs);
                foreach (DOFType dofType in element.ElementType.DOFEnumerator.GetDOFTypes(element)[i])
                {
                    bool isConstrainedDof = constrainedDOFs.TryGetValue(dofType, out double constraintDisplacement);
                    if (isConstrainedNode && isConstrainedDof)
                    {
                        Debug.Assert(elementNodalDisplacements[pos] == 0);
                        elementNodalDisplacements[pos] = constraintDisplacement;
                    }
                    pos++;
                }
            }
        }

        ////TODO: this should return Vector
        //public double[] GetLocalVectorFromGlobal(Element element, IVectorView globalVector) //TODOMaria: here is where the element displacements are assigned to zero if they are restrained
        //{                                                                                   //TODOMaria: Change visibility to private
        //    int localDOFs = 0;
        //    foreach (IList<DOFType> dofs in element.ElementType.DOFEnumerator.GetDOFTypes(element)) localDOFs += dofs.Count;
        //    var localVector = new double[localDOFs]; //TODOMaria: here is where I have to check if the dof is constrained

        //    int pos = 0;
        //    IList<IList<DOFType>> nodalDofs = element.ElementType.DOFEnumerator.GetDOFTypes(element);
        //    IList<INode> nodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
        //    for (int i = 0; i < nodes.Count; i++)
        //    {
        //        foreach (DOFType dofType in nodalDofs[i])
        //        {
        //            int dof = NodalDOFsDictionary[nodes[i].ID][dofType];
        //            if (dof != -1) localVector[pos] = globalVector[dof];
        //            pos++;
        //        }
        //    }
        //    return localVector;
        //}

        //public void AddLocalVectorToGlobal(Element element, double[] localVector, double[] globalVector)
        //{
        //    int pos = 0;
        //    IList<IList<DOFType>> nodalDofs = element.ElementType.DOFEnumerator.GetDOFTypes(element);
        //    IList<INode> nodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
        //    for (int i = 0; i < nodes.Count; i++)
        //    {
        //        foreach (DOFType dofType in nodalDofs[i])
        //        {
        //            int dof = NodalDOFsDictionary[nodes[i].ID][dofType];
        //            if (dof != -1) globalVector[dof] += localVector[pos];
        //            pos++;
        //        }
        //    }
        //}

        public IVector GetRHSFromSolution(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(DofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in elementsDictionary.Values)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria
                double[] localSolution = CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = CalculateElementNodalDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.MaterialsModified = true;
                var f = Vector.CreateFromArray(element.ElementType.CalculateForces(element, localSolution, localdSolution));
                DofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void SaveMaterialState()
        {
            foreach (Element element in elementsDictionary.Values) element.ElementType.SaveMaterialState();
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.MaterialsModified = false;
            foreach (Element element in elementsDictionary.Values) element.ElementType.ResetMaterialModified();
        }

        public void ClearMaterialStresses()
        {
            foreach (Element element in elementsDictionary.Values) element.ElementType.ClearMaterialStresses();
        }
    }
}
