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
        private readonly List<XNode2D> nodes;
        private readonly List<Element2D> elements;
        private readonly List<IEnrichmentItem2D> enrichments;

        // TODO: There should probably be a dedicated Constraint or Constraints or BoundaryCondition(s) class
        private readonly Dictionary<Node2D, SortedSet<StandardDOFType>> constraints;

        /// <summary>
        /// Multiple loads for the same dof are allowed for now. Probably they shouldn't.
        /// </summary>
        private readonly List<NodalLoad2D> loads;

        public IReadOnlyList<XNode2D> Nodes { get { return nodes; } }
        public IReadOnlyList<Element2D> Elements { get { return elements; } }
        public IReadOnlyList<IEnrichmentItem2D> Enrichments { get { return enrichments; } }
        public DOFEnumerator DofEnumerator { get; private set; }

        public Model2D()
        {
            this.nodes = new List<XNode2D>();
            this.elements = new List<Element2D>();
            this.enrichments = new List<IEnrichmentItem2D>();
            this.constraints = new Dictionary<Node2D, SortedSet<StandardDOFType>>();
            this.loads = new List<NodalLoad2D>();
        }

        public void AddNode(XNode2D node)
        {
            if (nodes.Contains(node))
                throw new ArgumentException("There is already a node with id = " + node.ID);
            nodes.Add(node);
        }

        public void AddElement(Element2D element)
        {
            if (elements.Contains(element))
                throw new ArgumentException("There is already a element with id = " + element.ID);
            elements.Add(element);
        }

        public void AddEnrichment(IEnrichmentItem2D enrichment)
        {
            if (enrichments.Contains(enrichment))
                throw new ArgumentException("There is already a enrichment with id = " + enrichment.ID);
            enrichments.Add(enrichment);
        }

        //TODO: Should I use the node's id instead? In a UI class, I probably should.
        public void AddConstraint(Node2D node, StandardDOFType dofType)
        {
            if (!nodes.Contains(node)) // TODO: This should be done more efficiently than O(N)
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
            if (!nodes.Contains(node)) // TODO: This should be done more efficiently than O(N)
            {
                throw new ArgumentException("There is no such node");
            }

            this.loads.Add(new NodalLoad2D(node, dofType, value));
        }

        public void EnumerateDofs()
        {
            this.DofEnumerator = new DOFEnumerator(nodes, constraints, elements);
        }

        public double[] CalculateForces()
        {
            if (DofEnumerator == null) EnumerateDofs();
            double[] rhs = new double[DofEnumerator.FreeStandardDofsCount + DofEnumerator.ArtificialDofsCount];
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
