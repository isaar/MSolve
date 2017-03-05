using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class Model2D
    {
        private Dictionary<int, XNode2D> nodes;
        private Dictionary<int, Element2D> elements;
        private Dictionary<int, IEnrichmentItem2D> enrichments;

        public DOFEnumerator DofEnumerator { get; private set; }

        public Model2D()
        {
            this.nodes = new Dictionary<int, XNode2D>();
            this.elements = new Dictionary<int, Element2D>();
            this.enrichments = new Dictionary<int, IEnrichmentItem2D>();
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

        public Dictionary<Node2D, StandardDOFType[]> FindConstraints()
        {
            // TODO: Implement or replace this method.
            return new Dictionary<Node2D, StandardDOFType[]>();
        }

        public void EnumerateDofs()
        {
            this.DofEnumerator = new DOFEnumerator(nodes.Values, FindConstraints(), elements.Values);
        }
    }
}
