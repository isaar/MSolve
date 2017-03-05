using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    class DOFEnumerator
    {
        private readonly Dictionary<Node2D, Dictionary<StandardDOFType, int>> standardDofs;
        private readonly Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>> artificialDofs;

        public DOFEnumerator(IEnumerable<XNode2D> nodes, Dictionary<Node2D, StandardDOFType[]> constraints,
            IEnumerable<Element2D> elements)
        {
            standardDofs = EnumerateStandardDofs(elements, constraints);
            artificialDofs = EnumerateArtificialDofs(nodes);
        }

        public IEnumerable<int> GetStandardDofSOfNode(Node2D node)
        {
            return standardDofs[node].Values;
        }

        public IEnumerable<int> GetArtificialDofSOfNode(XNode2D node)
        {
            return artificialDofs[node].Values;
        }

        private static Dictionary<Node2D, Dictionary<StandardDOFType, int>> EnumerateStandardDofs(
            IEnumerable<Element2D> elements, Dictionary<Node2D, StandardDOFType[]> constraints)
        {
            Dictionary<Node2D, HashSet<StandardDOFType>> nodalDOFTypes = FindUniqueDOFTypes(elements);

            var totalDofs = new Dictionary<Node2D, Dictionary<StandardDOFType, int>>();
            int dofCounter = 0;
            foreach (var pair in nodalDOFTypes)
            {
                Node2D node = pair.Key;
                var dofsOfThisNode = new Dictionary<StandardDOFType, int>();
                foreach (StandardDOFType dofType in pair.Value)
                {
                    if (constraints[node].Contains(dofType)) dofsOfThisNode[dofType] = -1; // Is this the best way? How about using the Ksf and Kss to calculate the reactions?
                    else dofsOfThisNode[dofType] = dofCounter++;
                }
                totalDofs[node] = dofsOfThisNode;
            }
            return totalDofs;
        }

        private static Dictionary<Node2D, HashSet<StandardDOFType>> FindUniqueDOFTypes(IEnumerable<Element2D> elements)
        {
            var totalDofs = new Dictionary<Node2D, HashSet<StandardDOFType>>();
            foreach (Element2D element in elements)
            {
                foreach (var entry in element.ElementType.StandardFiniteElement.GetNodalDOFTypes())
                {
                    HashSet<StandardDOFType> dofsOfThisNode;
                    bool alreadyExists = totalDofs.TryGetValue(entry.Key, out dofsOfThisNode);
                    if (!alreadyExists)
                    {
                        dofsOfThisNode = new HashSet<StandardDOFType>();
                        totalDofs.Add(entry.Key, dofsOfThisNode);
                    }
                    dofsOfThisNode.UnionWith(entry.Value);
                }
            }
            return totalDofs;
        }

        // Each artificial dof has index that is node major, then enrichment item major, then enrichment function major and finally axis minor
        private Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>> EnumerateArtificialDofs(IEnumerable<XNode2D> nodes)
        {
            var totalDofs = new Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>>();
            int dofCounter = standardDofs.Count; // This if I put everything in the same matrix
            //int dofCounter = 0; // This if I use different matrices
            foreach (XNode2D node in nodes)
            {
                var dofsOfThisNode = new Dictionary<ArtificialDOFType, int>();
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        dofsOfThisNode[dofType] = dofCounter++;
                    }
                }
                totalDofs[node] = dofsOfThisNode;
            }
            return totalDofs;
        }
    }
}
