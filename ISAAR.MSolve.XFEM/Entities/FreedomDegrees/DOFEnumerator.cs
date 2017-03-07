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

        public int StandardDofsCount { get { return standardDofs.Count; } }
        public int ArtificialDofsCount { get { return artificialDofs.Count; } }

        public DOFEnumerator(IEnumerable<XNode2D> nodes, Dictionary<Node2D, StandardDOFType[]> constraints,
            IEnumerable<Element2D> elements)
        {
            standardDofs = EnumerateStandardDofs(elements, constraints);
            artificialDofs = EnumerateArtificialDofs(nodes, standardDofs.Count);
        }

        public IEnumerable<int> GetStandardDofsOf(Node2D node)
        {
            return standardDofs[node].Values;
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Read is supposed to be O(1), but still...
        public int GetStandardDofOf(Node2D node, StandardDOFType dofType)
        {
            return standardDofs[node][dofType];
        }

        public IEnumerable<int> GetArtificialDofsOf(XNode2D node)
        {
            return artificialDofs[node].Values;
        }

        public int GetArtificialDofOf(XNode2D node, ArtificialDOFType dofType)
        {
            return artificialDofs[node][dofType];
        }

        public IReadOnlyList<int> MatchElementToGlobalStandardDofsOf(Element2D element)
        {
            var nodalDofTypes = element.ElementType.StandardFiniteElement.GetNodalDOFTypes();
            int[] globalDofs = new int[nodalDofTypes.Count];
            int elementDof = 0;
            foreach (KeyValuePair<Node2D, HashSet<StandardDOFType>> pair in nodalDofTypes)
            {
                Dictionary<StandardDOFType, int> globalDofsOfThisNode = standardDofs[pair.Key];
                // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                foreach (StandardDOFType dofType in pair.Value)
                {
                    globalDofs[elementDof] = globalDofsOfThisNode[dofType];
                    ++elementDof;
                }
            }
            return globalDofs;
        }

        public IReadOnlyList<int> MatchElementToGlobalArtificialDofsOf(Element2D element)
        {
            var globalDofs = new List<int>();
            //var nodalDofTypes = element.ElementType.StandardFiniteElement.GetNodalDOFTypes();
            
            int elementDof = 0;
            foreach (XNode2D node in element.ElementType.Nodes)
            {
                Dictionary<ArtificialDOFType, int> globalDofsOfThisNode = artificialDofs[node];
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) and drawn from XElement2D instead form enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        globalDofs[elementDof] = globalDofsOfThisNode[dofType];
                        ++elementDof;
                    }
                }
            }
            return globalDofs;
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
                    else dofsOfThisNode[dofType] = ++dofCounter;
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
        private Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>> EnumerateArtificialDofs(
            IEnumerable<XNode2D> nodes, int standardDofsCount)
        {
            var totalDofs = new Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>>();
            int dofCounter = standardDofsCount; // This if I put everything in the same matrix
            //int dofCounter = 0; // This if I use different matrices
            foreach (XNode2D node in nodes)
            {
                var dofsOfThisNode = new Dictionary<ArtificialDOFType, int>();
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        dofsOfThisNode[dofType] = ++dofCounter;
                    }
                }
                totalDofs[node] = dofsOfThisNode;
            }
            return totalDofs;
        }
    }
}
