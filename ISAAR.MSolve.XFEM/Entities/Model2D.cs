using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// This class both manages and assembles the FEM entities. TODO: Those 2 should be split. I could have a nested 
    /// Builder or better yet a UI class for assembling the model.
    /// </summary>
    class Model2D
    {
        private readonly Dictionary<int, XNode2D> nodes;
        private readonly Dictionary<int, Element2D> elements;
        private readonly Dictionary<int, IEnrichmentItem2D> enrichments;

        // TODO: There should probably be a dedicated Constraint or Constraints or BoundaryCondition(s) class
        private readonly Dictionary<Node2D, SortedSet<StandardDOFType>> constraints;

        /// <summary>
        /// Multiple loads for the same dof are allowed for now. Probably they shouldn't.
        /// </summary>
        private readonly List<NodalLoad2D> loads;

        public IEnumerable<XNode2D> Nodes { get { return nodes.Values; } }
        public IEnumerable<Element2D> Elements { get { return elements.Values; } }
        public IEnumerable<IEnrichmentItem2D> Enrichments { get { return enrichments.Values; } }
        public DOFEnumerator DofEnumerator { get; private set; }

        public Model2D()
        {
            this.nodes = new Dictionary<int, XNode2D>();
            this.elements = new Dictionary<int, Element2D>();
            this.enrichments = new Dictionary<int, IEnrichmentItem2D>();
            this.constraints = new Dictionary<Node2D, SortedSet<StandardDOFType>>();
            this.loads = new List<NodalLoad2D>();
        }

        public void AddNode(XNode2D node)
        {
            if (nodes.ContainsKey(node.ID))
                throw new ArgumentException("There is already a node with id = " + node.ID);
            nodes[node.ID] = node;
        }

        public void AddElement(Element2D element)
        {
            if (elements.ContainsKey(element.ID))
                throw new ArgumentException("There is already a element with id = " + element.ID);
            elements[element.ID] = element;
        }

        public void AddEnrichment(IEnrichmentItem2D enrichment)
        {
            if (enrichments.ContainsKey(enrichment.ID))
                throw new ArgumentException("There is already a enrichment with id = " + enrichment.ID);
            enrichments[enrichment.ID] = enrichment;
        }

        //TODO: Should I use the node's id instead? In a UI class, I probably should.
        public void AddConstraint(Node2D node, StandardDOFType dofType)
        {
            if (!nodes.Values.Contains(node)) // TODO: This should be done more efficiently than O(N)
            {
                throw new ArgumentException("There is no such node");
            }
             
            SortedSet<StandardDOFType> constraintsOfNode;
            bool alreadyExists = constraints.TryGetValue(node, out constraintsOfNode);
            if (!alreadyExists)
            {
                constraintsOfNode = new SortedSet<StandardDOFType>();
                constraints.Add(node, constraintsOfNode);
            }
            constraintsOfNode.Add(dofType);
        }

        //TODO: Should I use the node's id instead? In a UI class, I probably should.
        public void AddNodalLoad(Node2D node, StandardDOFType dofType, double value)
        {
            if (!nodes.Values.Contains(node)) // TODO: This should be done more efficiently than O(N)
            {
                throw new ArgumentException("There is no such node");
            }

            this.loads.Add(new NodalLoad2D(node, dofType, value));
        }

        public void EnumerateDofs()
        {
            this.DofEnumerator = new DOFEnumerator(nodes.Values, constraints, elements.Values);
        }

        public double[] CalculateForces()
        {
            if (DofEnumerator == null) EnumerateDofs();
            double[] rhs = new double[DofEnumerator.FreeStandardDofsCount];
            foreach (NodalLoad2D load in loads)
            {
                int dof = DofEnumerator.GetStandardDofOf(load.Node, load.DofType);
                if (dof < 0) throw new NotImplementedException("Load on a constraint dof at node "
                    + load.Node.ID + ", axis " + load.DofType);
                else rhs[dof] += load.Value;
            }
            return rhs;
        }
    }
}
