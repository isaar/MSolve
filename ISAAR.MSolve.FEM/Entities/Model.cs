using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Model:IStructuralModel
    {
        public const int constrainedDofIdx = -1;
        private int totalDOFs = 0;
        private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        //private readonly IList<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Subdomain> subdomainsDictionary = new Dictionary<int, Subdomain>();
        private readonly Dictionary<int, Cluster> clustersDictionary = new Dictionary<int, Cluster>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly IList<Load> loads = new List<Load>();
        private readonly IList<ElementMassAccelerationLoad> elementMassAccelerationLoads = new List<ElementMassAccelerationLoad>();
        private readonly IList<MassAccelerationLoad> massAccelerationLoads = new List<MassAccelerationLoad>();
        private readonly IList<IMassAccelerationHistoryLoad> massAccelerationHistoryLoads = new List<IMassAccelerationHistoryLoad>();
        private readonly IList<ElementMassAccelerationHistoryLoad> elementMassAccelerationHistoryLoads = new List<ElementMassAccelerationHistoryLoad>();

        #region Properties
        //public IList<EmbeddedNode> EmbeddedNodes
        //{
        //    get { return embeddedNodes; }
        //}

        public Dictionary<int, Node> NodesDictionary
        {
            get { return nodesDictionary; }
        }

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, Subdomain> SubdomainsDictionary
        {
            get { return subdomainsDictionary; }
        }

	    public Dictionary<int, ISubdomain> ISubdomainsDictionary
	    {
		    get
		    {
			    var a = new Dictionary<int, ISubdomain>();
			    foreach (var subdomain in subdomainsDictionary.Values)
				    a.Add(subdomain.ID,subdomain);
				return a;
		    }
	    }

		public Dictionary<int, Cluster> ClustersDictionary
        {
            get { return clustersDictionary; }
        }

        public IList<Node> Nodes
        {
            get { return nodesDictionary.Values.ToList<Node>(); }
        }

        public IList<Element> Elements
        {
            get { return elementsDictionary.Values.ToList<Element>(); }
        }

        public IList<Subdomain> Subdomains
        {
            get { return subdomainsDictionary.Values.ToList<Subdomain>(); }
        }

        public IList<Cluster> Clusters
        {
            get { return clustersDictionary.Values.ToList<Cluster>(); }
        }

        public IList<Load> Loads
        {
            get { return loads; }
        }

        public IList<ElementMassAccelerationLoad> ElementMassAccelerationLoads
        {
            get { return elementMassAccelerationLoads; }
        }

        public IList<MassAccelerationLoad> MassAccelerationLoads
        {
            get { return massAccelerationLoads; }
        }

        public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads
        {
            get { return massAccelerationHistoryLoads; }
        }

        public IList<ElementMassAccelerationHistoryLoad> ElementMassAccelerationHistoryLoads
        {
            get { return elementMassAccelerationHistoryLoads; }
        }

        public Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary
        {
            get { return nodalDOFsDictionary; }
        }

        public int TotalDOFs
        {
            get { return totalDOFs; }
        }

        #endregion

        #region Data interconnection routines
        private void BuildElementDictionaryOfEachNode()
        {
            foreach (Element element in elementsDictionary.Values)
                foreach (Node node in element.Nodes)
                    node.ElementsDictionary.Add(element.ID, element);
        }

        private void BuildSubdomainOfEachElement()
        {
            foreach (Subdomain subdomain in subdomainsDictionary.Values)
                foreach (Element element in subdomain.ElementsDictionary.Values)
                    element.Subdomain = subdomain;
        }

        private void BuildInterconnectionData()
        {
            BuildSubdomainOfEachElement();
            DuplicateInterSubdomainEmbeddedElements();
            BuildElementDictionaryOfEachNode();
            foreach (Node node in nodesDictionary.Values)
                node.BuildSubdomainDictionary();

            //BuildNonConformingNodes();

            foreach (Subdomain subdomain in subdomainsDictionary.Values)
                subdomain.BuildNodesDictionary();
        }

        private void DuplicateInterSubdomainEmbeddedElements()
        {
            foreach (var e in elementsDictionary.Values.Where(x => x.ElementType is IEmbeddedElement))
            {
                var subs = ((IEmbeddedElement)e.ElementType).EmbeddedNodes.Select(x => x.EmbeddedInElement.Subdomain).Distinct();
                foreach (var s in subs.Where(x => x.ID != e.Subdomain.ID))
                    s.ElementsDictionary.Add(e.ID, e);
            }
        }

        private void BuildNonConformingNodes()
        {
            List<int> subIDs = new List<int>();
            foreach (Element element in elementsDictionary.Values)
            {
                subIDs.Clear();

                foreach (Node node in element.Nodes)
                    foreach (int subID in node.SubdomainsDictionary.Keys)
                        if (!subIDs.Contains(subID)) subIDs.Add(subID);

                foreach (Node node in element.Nodes)
                    foreach (int subID in subIDs)
                        if (!node.SubdomainsDictionary.ContainsKey(subID)) node.NonMatchingSubdomainsDictionary.Add(subID, subdomainsDictionary[subID]);
            }
        }

        private void EnumerateGlobalDOFs()
        {
            totalDOFs = 0;
            Dictionary<int, List<DOFType>> nodalDOFTypesDictionary = new Dictionary<int, List<DOFType>>();
            foreach (Element element in elementsDictionary.Values)
            {
                for (int i = 0; i < element.Nodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(element.Nodes[i].ID))
                        nodalDOFTypesDictionary.Add(element.Nodes[i].ID, new List<DOFType>());
                    nodalDOFTypesDictionary[element.Nodes[i].ID].AddRange(element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            foreach (Node node in nodesDictionary.Values)
            {
                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                foreach (DOFType dofType in nodalDOFTypesDictionary[node.ID].Distinct<DOFType>())
                {
                    int dofID = 0;
                    foreach (DOFType constraint in node.Constraints)
                    {
                        if (constraint == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }

                    //// TODO: this is not applicable! Embedded nodes have to do with the embedded element and not with the host
                    //// User should define which DOFs would be dependent on the host element. For our case
                    //// we should select between translational and translational + rotational
                    //var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                    ////if (node.EmbeddedInElement != null && node.EmbeddedInElement.ElementType.GetDOFTypes(null)
                    ////    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    ////    dofID = -1;
                    //if (embeddedNode != null && embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null)
                    //    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    //    dofID = -1;

                    if (dofID == 0)
                    {
                        dofID = totalDOFs;
                        totalDOFs++;
                    }
                    dofsDictionary.Add(dofType, dofID);
                }

                nodalDOFsDictionary.Add(node.ID, dofsDictionary);
            }
        }

        private void EnumerateDOFs()
        {
            EnumerateGlobalDOFs();
            foreach (Subdomain subdomain in subdomainsDictionary.Values)
            {
                subdomain.EnumerateDOFs();
                subdomain.AssignGlobalNodalDOFsFromModel(nodalDOFsDictionary);
            }
        }

        private void AssignNodalLoads()
        {
            foreach (Subdomain subdomain in subdomainsDictionary.Values)
                Array.Clear(subdomain.Forces, 0, subdomain.Forces.Length);
            foreach (Load load in loads)
                foreach (Subdomain subdomain in load.Node.SubdomainsDictionary.Values)
                {
                    int dof = subdomain.NodalDOFsDictionary[load.Node.ID][load.DOF];
                    if (dof >= 0)
                        subdomain.Forces[dof] = load.Amount / load.Node.SubdomainsDictionary.Count;
                }
        }

        private void AssignElementMassLoads()
        {
            foreach (ElementMassAccelerationLoad load in elementMassAccelerationLoads)
                load.Element.Subdomain.AddLocalVectorToGlobal(load.Element,
                    load.Element.ElementType.CalculateAccelerationForces(load.Element, massAccelerationLoads),
                    load.Element.Subdomain.Forces);
        }

        private void AssignMassAccelerationLoads()
        {
            if (massAccelerationLoads.Count < 1) return;

            foreach (Subdomain subdomain in subdomainsDictionary.Values)
                foreach (Element element in subdomain.ElementsDictionary.Values)
                    subdomain.AddLocalVectorToGlobal(element,
                        element.ElementType.CalculateAccelerationForces(element, massAccelerationLoads),
                        subdomain.Forces);
        }

        public void AssignLoads()
        {
            AssignNodalLoads();
            AssignElementMassLoads();
            AssignMassAccelerationLoads();
        }

        public void AssignMassAccelerationHistoryLoads(int timeStep)
        {
            if (massAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(massAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in massAccelerationHistoryLoads)
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });

                foreach (Subdomain subdomain in subdomainsDictionary.Values)
                    foreach (Element element in subdomain.ElementsDictionary.Values)
                        subdomain.AddLocalVectorToGlobal(element,
                            element.ElementType.CalculateAccelerationForces(element, m), subdomain.Forces);
            }

            foreach (ElementMassAccelerationHistoryLoad load in elementMassAccelerationHistoryLoads)
            {
                MassAccelerationLoad hl = new MassAccelerationLoad() { Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF };
                load.Element.Subdomain.AddLocalVectorToGlobal(load.Element,
                    load.Element.ElementType.CalculateAccelerationForces(load.Element, (new MassAccelerationLoad[] { hl }).ToList()),
                    load.Element.Subdomain.Forces);
            }
        }

        public void ConnectDataStructures()
        {
            BuildInterconnectionData();
            EnumerateDOFs();
            //EnumerateSubdomainLagranges();
            //EnumerateDOFMultiplicity();
            AssignLoads();
        }
        #endregion

        public void Clear()
        {
            loads.Clear();
            clustersDictionary.Clear();
            subdomainsDictionary.Clear();
            elementsDictionary.Clear();
            nodesDictionary.Clear();
            nodalDOFsDictionary.Clear();
            elementMassAccelerationHistoryLoads.Clear();
            elementMassAccelerationLoads.Clear();
            massAccelerationHistoryLoads.Clear();
            massAccelerationLoads.Clear();
        }
    }
}
