using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Enrichments.Items;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    /// TODO: abstract this with an IGlobalEnumerator2D interface. That way I could have different dof orders:
    /// e.g. 1) first all standard dofs and then all enriched, 2) First all dofs (std & enr) of the 1st node, 
    /// then of the 2nd, etc.
    /// TODO: there is no reason to interact with Element2D instead of directly XContinuumElement. I do not use 
    /// anything from Element2D and accessing element.ElementType.SomeMethod() violates Diameter's law.
    class DOFEnumerator
    {
        private readonly int freeStandardDofsCount;
        private readonly int artificialDofsCount;

        private readonly Dictionary<Node2D, Dictionary<StandardDOFType, int>> standardDofs;
        private readonly Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>> artificialDofs;

        public int FreeStandardDofsCount { get { return freeStandardDofsCount; } }
        public int ConstrainedStandardDofsCount { get; }
        public int ArtificialDofsCount { get { return artificialDofsCount; } }
        

        // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections used to represent constraints
        public DOFEnumerator(IEnumerable<XNode2D> nodes, Dictionary<Node2D, SortedSet<StandardDOFType>> constraints,
            IEnumerable<Element2D> elements)
        {
            EnumerateStandardDofs(elements, constraints, out freeStandardDofsCount, out standardDofs);
            EnumerateArtificialDofs(nodes, freeStandardDofsCount, out artificialDofsCount, out artificialDofs);

            ConstrainedStandardDofsCount = 0;
            foreach (SortedSet<StandardDOFType> dofs in constraints.Values) ConstrainedStandardDofsCount += dofs.Count;
        }

        public IEnumerable<int> GetStandardDofsOf(Node2D node)
        {
            return standardDofs[node].Values;
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Dictionary's Read operation is supposed to be O(1), but still...
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
            var nodalDofTypes = element.ElementType.GetStandardNodalDOFTypes();
            var globalDofs = new List<int>(); // An array might be quicker, but it requires counting the elements' dofs
            foreach (KeyValuePair<Node2D, HashSet<StandardDOFType>> pair in nodalDofTypes)
            {
                Dictionary<StandardDOFType, int> globalDofsOfThisNode = standardDofs[pair.Key];
                // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                foreach (StandardDOFType dofType in pair.Value)
                {
                    globalDofs.Add(globalDofsOfThisNode[dofType]);
                }
            }
            return globalDofs;
        }

        public IReadOnlyList<int> MatchElementToGlobalArtificialDofsOf(Element2D element)
        {
            var globalDofs = new List<int>();
            foreach (XNode2D node in element.ElementType.Nodes)
            {
                Dictionary<ArtificialDOFType, int> globalDofsOfThisNode = artificialDofs[node];
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) and drawn from XElement2D instead form enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        globalDofs.Add(globalDofsOfThisNode[dofType]);
                    }
                }
            }
            return globalDofs;
        }

        /// <summary>
        /// For constrained dofs the corresponding displacement is 0.0
        /// TODO: Modify this method to extract any kind of vector, not only displacement, which means different 
        /// handling of constrained dofs, if constrained dofs are defined in the first place.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="globalVector"></param>
        /// <returns></returns>
        public double[] ExtractDisplacementVectorOfElementFromGlobal(Element2D element, double[] globalVector)
        {
            IReadOnlyDictionary<Node2D, HashSet<StandardDOFType>> elementDofs = 
                element.ElementType.GetStandardNodalDOFTypes();
            int dofsCount = 0;
            foreach (HashSet<StandardDOFType> nodalDofs in elementDofs.Values) dofsCount += nodalDofs.Count;
            double[] localVector = new double[dofsCount];

            int localDof = 0;
            foreach (KeyValuePair<Node2D, HashSet<StandardDOFType>> pair in elementDofs)
            {
                Dictionary<StandardDOFType, int> globalDofTypesOfNode = standardDofs[pair.Key];
                foreach (StandardDOFType localDofType in pair.Value)
                {
                    int globalDof = globalDofTypesOfNode[localDofType];
                    if (globalDof == -1) localVector[localDof++] = 0;
                    else localVector[localDof++] = globalVector[globalDof];
                }
            }
            return localVector;
        }

        private static void EnumerateStandardDofs(
            IEnumerable<Element2D> elements, Dictionary<Node2D, SortedSet<StandardDOFType>> constraints, 
            out int standardDofsCount, out Dictionary<Node2D, Dictionary<StandardDOFType, int>> standardDofs)
        {
            IDictionary<Node2D, HashSet<StandardDOFType>> nodalDOFTypes = FindUniqueDOFTypes(elements);

            standardDofs = new Dictionary<Node2D, Dictionary<StandardDOFType, int>>();
            int dofCounter = 0;
            foreach (var pair in nodalDOFTypes)
            {
                Node2D node = pair.Key;

                SortedSet<StandardDOFType> constrainedDofs;
                bool hasConstraints = constraints.TryGetValue(node, out constrainedDofs);
                if (!hasConstraints) constrainedDofs = new SortedSet<StandardDOFType>(); // empty constraints

                var dofsOfThisNode = new Dictionary<StandardDOFType, int>();
                foreach (StandardDOFType dofType in pair.Value)
                {
                    // Perhaps I could have a function that splits the dofs into free/constrained, instead of looking up each dof in the constrained ones.
                    if (constrainedDofs.Contains(dofType)) dofsOfThisNode[dofType] = -1; // Is this the best way? How about using the Ksf and Kss to calculate the reactions?
                    else dofsOfThisNode[dofType] = dofCounter++;
                }
                standardDofs[node] = dofsOfThisNode;
            }
            standardDofsCount = dofCounter;
        }

        private static IDictionary<Node2D, HashSet<StandardDOFType>> FindUniqueDOFTypes(IEnumerable<Element2D> elements)
        {
            var totalDofs = new SortedDictionary<Node2D, HashSet<StandardDOFType>>();
            foreach (Element2D element in elements)
            {
                foreach (var entry in element.ElementType.GetStandardNodalDOFTypes())
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
        private void EnumerateArtificialDofs(IEnumerable<XNode2D> nodes, int standardDofsCount,
            out int artificialDofsCount, out Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>> artificialDofs)
        {
            artificialDofs = new Dictionary<XNode2D, Dictionary<ArtificialDOFType, int>>();
            int dofCounter = standardDofsCount; // This if I put everything in the same matrix
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
                artificialDofs[node] = dofsOfThisNode;
            }
            artificialDofsCount = dofCounter - standardDofsCount;
        }
    }
}
