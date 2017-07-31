using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// This class both manages and assembles the FEM entities. TODO: Those 2 should be split. I could have a nested 
    /// Builder or better yet a UI class for assembling the model.
    /// </summary>
    class Model2D
    {
        private readonly List<XNode2D> nodes;
        private readonly List<XContinuumElement2D> elements;
        private readonly List<IEnrichmentItem2D> enrichments;
        private readonly Table<XNode2D, StandardDOFType, double> constraints;

        /// <summary>
        /// Multiple loads for the same dof are not allowed. Any attempt to input them will result in an exception 
        /// thrown by the table object.
        /// </summary>
        private readonly Table<XNode2D, StandardDOFType, double> loads;

        public IReadOnlyList<XNode2D> Nodes { get { return nodes; } }
        public IReadOnlyList<XContinuumElement2D> Elements { get { return elements; } }
        public IReadOnlyList<IEnrichmentItem2D> Enrichments { get { return enrichments; } }
        public ITable<XNode2D, StandardDOFType, double> Constraints { get { return constraints; } }
        public DOFEnumerator DofEnumerator { get; private set; }

        public Model2D()
        {
            this.nodes = new List<XNode2D>();
            this.elements = new List<XContinuumElement2D>();
            this.enrichments = new List<IEnrichmentItem2D>();
            this.constraints = new Table<XNode2D, StandardDOFType, double>();
            this.loads = new Table<XNode2D, StandardDOFType, double>();
        }

        public void AddNode(XNode2D node)
        {
            if (nodes.Contains(node))
                throw new ArgumentException("There is already a node with id = " + node.ID);
            nodes.Add(node);
        }

        public void AddElement(XContinuumElement2D element)
        {
            if (elements.Contains(element))
                throw new ArgumentException("This element is already inserted");
            elements.Add(element);
        }

        public void AddEnrichment(IEnrichmentItem2D enrichment)
        {
            if (enrichments.Contains(enrichment))
                throw new ArgumentException("This enrichment is already inserted");
            enrichments.Add(enrichment);
        }

        public void AddConstraint(XNode2D node, StandardDOFType dof, double displacement)
        {
            // TODO: This should be done more efficiently than O(N)
            if (!nodes.Contains(node)) throw new ArgumentException("There is no such node");
            constraints[node, dof] = displacement;
        }

        //TODO: Should I use the node's id instead? In a UI class, I probably should.
        public void AddNodalLoad(XNode2D node, StandardDOFType dof, double magnitude)
        {
            // TODO: This should be done more efficiently than O(N)
            if (!nodes.Contains(node)) throw new ArgumentException("There is no such node");
            loads[node, dof] = magnitude;
        }

        public void EnumerateDofs()
        {
            this.DofEnumerator = new DOFEnumerator(nodes, constraints, elements);
        }

        public double[] CalculateFreeForces()
        {
            if (DofEnumerator == null) EnumerateDofs();
            double[] rhs = new double[DofEnumerator.FreeDofsCount + DofEnumerator.ArtificialDofsCount];
            foreach (Tuple<XNode2D, StandardDOFType, double> entry in loads)
            {
                try
                {
                    int dof = DofEnumerator.GetFreeDofOf(entry.Item1, entry.Item2);
                    rhs[dof] += entry.Item3; // This supports multiple loads on the same dof, which that isn't implemented yet
                }
                catch (KeyNotFoundException ex)
                {
                    throw new NotImplementedException("Load on a constrained dof at node "
                    + entry.Item1.ID + ", axis " + entry.Item2, ex);
                }
            }
            return rhs;
        }

        public double[] CalculateConstrainedDisplacements()
        {
            if (DofEnumerator == null) EnumerateDofs();
            double[] uc = new double[DofEnumerator.ConstrainedDofsCount];
            foreach (Tuple<XNode2D, StandardDOFType, double> entry in constraints)
            {
                int dof = DofEnumerator.GetConstrainedDofOf(entry.Item1, entry.Item2);
                uc[dof] = entry.Item3;
            }
            return uc;
        }
    }
}
