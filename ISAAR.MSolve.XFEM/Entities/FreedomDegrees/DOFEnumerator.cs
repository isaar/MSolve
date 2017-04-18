using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    /// TODO: abstract this with an IGlobalEnumerator2D interface. That way I could have different dof orders:
    /// e.g. 1) first all standard dofs and then all enriched, 2) First all dofs (std & enr) of the 1st node, 
    /// then of the 2nd, etc.
    class DOFEnumerator
    {
        private readonly int freeStandardDofsCount;
        private readonly int artificialDofsCount;
        
        private readonly Table<XNode2D, StandardDOFType, int> standardDofs;

        // Nodes that have no artificial dofs, are not stored
        private readonly Table<XNode2D, ArtificialDOFType, int> artificialDofs;

        public int FreeStandardDofsCount { get { return freeStandardDofsCount; } }
        public int ConstrainedStandardDofsCount { get; }
        public int ArtificialDofsCount { get { return artificialDofsCount; } }
        

        // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections used to represent constraints
        public DOFEnumerator(IEnumerable<XNode2D> nodes, ITable<XNode2D, StandardDOFType, double> constraints,
            IEnumerable<XContinuumElement2D> elements)
        {
            EnumerateStandardDofs(elements, constraints, out freeStandardDofsCount, out standardDofs);
            EnumerateArtificialDofs(nodes, freeStandardDofsCount, out artificialDofsCount, out artificialDofs);

            ConstrainedStandardDofsCount = constraints.EntryCount;
        }

        public IEnumerable<int> GetStandardDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            return standardDofs.GetValuesOfRow(node);
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Dictionary's Read operation is supposed to be O(1), but still...
        public int GetStandardDofOf(XNode2D node, StandardDOFType dofType)
        {
            return standardDofs[node, dofType];
        }

        public IEnumerable<int> GetArtificialDofsOf(XNode2D node)
        {
            try
            {
                return artificialDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException ex) 
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public int GetArtificialDofOf(XNode2D node, ArtificialDOFType dofType)
        {
            return artificialDofs[node, dofType];
        }

        public IReadOnlyList<int> MatchElementToGlobalStandardDofsOf(XContinuumElement2D element)
        {
            var nodalDofTypes = element.GetStandardNodalDOFTypes();
            var globalDofs = new List<int>(); // An array might be quicker, but it requires counting the elements' dofs
            foreach (KeyValuePair<XNode2D, HashSet<StandardDOFType>> pair in nodalDofTypes)
            {
                // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                foreach (StandardDOFType dofType in pair.Value)
                {
                    globalDofs.Add(standardDofs[pair.Key, dofType]);
                }
            }
            return globalDofs;
        }

        public IReadOnlyList<int> MatchElementToGlobalArtificialDofsOf(XContinuumElement2D element)
        {
            var globalDofs = new List<int>();
            foreach (XNode2D node in element.Nodes)
            {
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) and drawn from XXContinuumElement2D instead form enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        globalDofs.Add(artificialDofs[node, dofType]);
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
        public double[] ExtractDisplacementVectorOfElementFromGlobal(XContinuumElement2D element, double[] globalVector)
        {
            IReadOnlyDictionary<XNode2D, HashSet<StandardDOFType>> elementDofs = 
                element.GetStandardNodalDOFTypes();
            int dofsCount = 0;
            foreach (HashSet<StandardDOFType> nodalDofs in elementDofs.Values) dofsCount += nodalDofs.Count;
            double[] localVector = new double[dofsCount];

            int localDof = 0;
            foreach (KeyValuePair<XNode2D, HashSet<StandardDOFType>> pair in elementDofs)
            {
                foreach (StandardDOFType localDofType in pair.Value)
                {
                    int globalDof = standardDofs[pair.Key, localDofType];
                    if (globalDof == -1) localVector[localDof++] = 0;
                    else localVector[localDof++] = globalVector[globalDof];
                }
            }
            return localVector;
        }

        private static void EnumerateStandardDofs(
            IEnumerable<XContinuumElement2D> elements, ITable<XNode2D, StandardDOFType, double> constraints, 
            out int standardDofsCount, out Table<XNode2D, StandardDOFType, int> standardDofs)
        {
            IDictionary<XNode2D, HashSet<StandardDOFType>> nodalDOFTypes = FindUniqueDOFTypes(elements);

            standardDofs = new Table<XNode2D, StandardDOFType, int>();
            int dofCounter = 0;
            foreach (var pair in nodalDOFTypes)
            {
                XNode2D node = pair.Key;
                foreach (StandardDOFType dofType in pair.Value)
                {
                    // Perhaps I could have a function that splits the dofs into free/constrained, instead of looking up each dof in the constrained ones.
                    if (constraints.Contains(node, dofType)) standardDofs[node, dofType] = -1; // Is this the best way? How about using the Ksf and Kss to calculate the reactions?
                    else standardDofs[node, dofType] = dofCounter++;
                }
            }
            standardDofsCount = dofCounter;
        }

        private static IDictionary<XNode2D, HashSet<StandardDOFType>> FindUniqueDOFTypes(IEnumerable<XContinuumElement2D> elements)
        {
            var totalDofs = new SortedDictionary<XNode2D, HashSet<StandardDOFType>>();
            foreach (XContinuumElement2D element in elements)
            {
                foreach (var entry in element.GetStandardNodalDOFTypes())
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
            out int artificialDofsCount, out Table<XNode2D, ArtificialDOFType, int> artificialDofs)
        {
            artificialDofs = new Table<XNode2D, ArtificialDOFType, int>();
            int dofCounter = standardDofsCount; // This if I put everything in the same matrix
            //int dofCounter = 0; // This if I use different matrices
            foreach (XNode2D node in nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        artificialDofs[node, dofType] = dofCounter++;
                    }
                }
            }
            artificialDofsCount = dofCounter - standardDofsCount;
        }
    }
}
