using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;

//TODO: There is a lot of repetition between this FEM.Model and IGA.Model with regards to interconnection data. That code should 
//      be moved to a common class. Same goes for the interconnection methods of XSubdomain.
namespace ISAAR.MSolve.XFEM.Entities
{
    public class XModel : IStructuralModel
    {
        public IDomain2DBoundary Boundary { get; set; }

        public Table<INode, IDofType, double> Constraints { get; private set; } = new Table<INode, IDofType, double>();

        IReadOnlyList<IElement> IStructuralModel.Elements => Elements;
        public List<IXFiniteElement> Elements { get; } = new List<IXFiniteElement>();

        public IGlobalFreeDofOrdering GlobalDofOrdering { get; set; }

        public List<NodalLoad> NodalLoads { get; private set; } = new List<NodalLoad>();

        public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads => throw new NotImplementedException();

        IReadOnlyList<INode> IStructuralModel.Nodes => Nodes;
        public List<XNode> Nodes { get; } = new List<XNode>();

        IReadOnlyList<ISubdomain> IStructuralModel.Subdomains => Subdomains.Values.ToList();
        public Dictionary<int, XSubdomain> Subdomains { get; } = new Dictionary<int, XSubdomain>();

        public void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
        {
            foreach (XSubdomain subdomain in Subdomains.Values) subdomain.Forces.Clear();
            AssignNodalLoads(distributeNodalLoads);
        }

        public void AssignNodalLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
        {
            var globalNodalLoads = new Table<INode, IDofType, double>();
            foreach (NodalLoad load in NodalLoads) globalNodalLoads.TryAdd(load.Node, load.DofType, load.Value);

            Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
            foreach (var idSubdomainLoads in subdomainNodalLoads)
            {
                Subdomains[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
            }
        }

        public void AssignMassAccelerationHistoryLoads(int timeStep) => throw new NotImplementedException();

        public void ConnectDataStructures()
        {
            BuildInterconnectionData();
            AssignConstraints();
            RemoveInactiveNodalLoads();
        }

        private void AssignConstraints()
        {
            foreach (XNode node in Nodes)
            {
                if (node.Constraints == null) continue;
                foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            }

            foreach (XSubdomain subdomain in Subdomains.Values) subdomain.ExtractConstraintsFromGlobal(Constraints);
        }

        private void BuildInterconnectionData()
        {
            BuildSubdomainOfEachElement();

            // Storing the elements of each node is done by the IMesh class, if necessary. TODO: Find out what problems this causes.
            BuildElementDictionaryOfEachNode();

            // TODO: Storing the subdomains of each node should be done by another class, if necessary.
            foreach (XNode node in Nodes) node.BuildXSubdomainDictionary();

            foreach (XSubdomain subdomain in Subdomains.Values) subdomain.DefineNodesFromElements();
        }

        private void BuildElementDictionaryOfEachNode()
        {
            foreach (IXFiniteElement element in Elements)
            {
                foreach (XNode node in element.Nodes) node.ElementsDictionary[element.ID] = element;
            }
        }

        private void BuildSubdomainOfEachElement()
        {
            foreach (XSubdomain subdomain in Subdomains.Values)
            {
                foreach (IXFiniteElement element in subdomain.Elements) element.Subdomain = subdomain;
            }
        }

        private void RemoveInactiveNodalLoads()
        {
            // Static loads
            var activeLoadsStatic = new List<NodalLoad>(NodalLoads.Count);
            foreach (NodalLoad load in NodalLoads)
            {
                bool isConstrained = Constraints.Contains(load.Node, load.DofType);
                if (!isConstrained) activeLoadsStatic.Add(load);
            }
            NodalLoads = activeLoadsStatic;
        }
    }
}
