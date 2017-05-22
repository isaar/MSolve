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
    /// TODO: The enumeration should be in a dedicated class DOFEnumerator. Managing and exposing the dof tables and 
    /// the conversion methods should be in a DOFManager class or in Model. The DOFManager will greatly violate LoD 
    /// (without wrapping access to the tables and most methods), while if they are in Model it will violate SRP
    /// TODO: The MatchElementToGlobal... methods do a lot of duplicate iterations over the dofs. Perhaps a completely different way 
    /// should be found, which also changes the way the GlobalAssembler works. The ideal would be to split each element 
    /// matrix into the free, constrained and artificial parts required by the Assembler, without any overlaps, which 
    /// are currently denoted as -1 and need to be checked every single time they are created and accessed.
    class DOFEnumerator
    {
        private readonly int freeDofsCount;
        private readonly int constrainedDofsCount;
        private readonly int artificialDofsCount;
        
        private readonly Table<XNode2D, StandardDOFType, int> freeDofs;
        private readonly Table<XNode2D, StandardDOFType, int> constrainedDofs;
        private readonly Table<XNode2D, ArtificialDOFType, int> artificialDofs;

        public int FreeDofsCount { get { return freeDofsCount; } }
        public int ConstrainedDofsCount { get { return constrainedDofsCount; } }
        public int ArtificialDofsCount { get { return artificialDofsCount; } }
        

        // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections used to represent constraints
        public DOFEnumerator(IEnumerable<XNode2D> nodes, ITable<XNode2D, StandardDOFType, double> constraints,
            IEnumerable<XContinuumElement2D> elements)
        {
            IDictionary<XNode2D, HashSet<StandardDOFType>> nodalDOFTypes = FindUniqueDOFTypes(elements);
            EnumerateFreeDofs(nodalDOFTypes, constraints, out freeDofsCount, out freeDofs);
            EnumerateConstrainedDofs(nodalDOFTypes, constraints, out constrainedDofsCount, out constrainedDofs);
            EnumerateArtificialDofs(nodes, freeDofsCount, out artificialDofsCount, out artificialDofs);
        }

        public IEnumerable<int> GetFreeDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            try
            {
                return freeDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException ex)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        // Would it be faster to return the Dictionary<StandardDofType, int> for consecutive accesses of the dofs of this node? 
        // Dictionary's Read operation is supposed to be O(1), but still...
        public int GetFreeDofOf(XNode2D node, StandardDOFType dofType)
        {
            return freeDofs[node, dofType];
        }

        public IEnumerable<int> GetConstrainedDofsOf(XNode2D node)
        {
            // Perhaps it would be more efficient for the client to traverse all (dofType, dofIdx) pairs 
            // than assembling the dofIdx collection beforehand.
            try
            {
                return constrainedDofs.GetValuesOfRow(node);
            }
            catch (KeyNotFoundException ex)
            {
                return new int[] { }; // TODO: It would be better if this method is not called for a non enriched node.
            }
        }

        public int GetConstrainedDofOf(XNode2D node, StandardDOFType dofType)
        {
            return constrainedDofs[node, dofType];
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

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof. 
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public void MatchElementToGlobalStandardDofsOf(XContinuumElement2D element, 
            out IReadOnlyDictionary<int, int> elementToGlobalFreeDofs, 
            out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs)
        {
            ITable<XNode2D, StandardDOFType, int> elementDofs = element.GetStandardDofs();
            var globalFreeDofs = new Dictionary<int, int>();
            var globalConstrainedDofs = new Dictionary<int, int>();
            foreach (Tuple<XNode2D, StandardDOFType, int> entry in elementDofs)
            {
                int freeGlobalDof;
                bool isFree = this.freeDofs.TryGetValue(entry.Item1, entry.Item2, out freeGlobalDof);
                if (isFree) globalFreeDofs[entry.Item3] = freeGlobalDof;
                else globalConstrainedDofs[entry.Item3] = this.constrainedDofs[entry.Item1, entry.Item2];
            }
            elementToGlobalFreeDofs = globalFreeDofs;
            elementToGlobalConstrainedDofs = globalConstrainedDofs;
        }

        /// <summary>
        /// Index i = element local dof. Dictionary[i] = global dof.
        /// </summary>
        /// <param name="element"></param>
        /// <returns></returns>
        public IReadOnlyDictionary<int, int> MatchElementToGlobalArtificialDofsOf(XContinuumElement2D element)
        {
            var globalArtificialDofs = new Dictionary<int, int>();
            int elementDof = 0;
            foreach (XNode2D node in element.Nodes)
            {
                // TODO: Perhaps the nodal dof types should be decided by the element type (structural, continuum) 
                // and drawn from XXContinuumElement2D instead of from enrichment items
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    // TODO: Perhaps I could iterate directly on the dofs, ignoring dof types for performance, if the order is guarranteed
                    foreach (ArtificialDOFType dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        globalArtificialDofs[elementDof++] = this.artificialDofs[node, dofType];
                    }
                }
            }
            return globalArtificialDofs;
        }

        /// <summary>
        /// TODO: Modify this method to extract any kind of vector, not only displacement, which means different 
        /// handling of constrained dofs, if constrained dofs are defined in the first place.
        /// </summary>
        /// <param name="element"></param>
        /// <param name="globalFreeVector"></param>
        /// <param name="globalConstrainedVector"></param>
        /// <returns></returns>
        public double[] ExtractDisplacementVectorOfElementFromGlobal(XContinuumElement2D element, 
            double[] globalFreeVector, double[] globalConstrainedVector)
        {
            ITable<XNode2D, StandardDOFType, int> elementDofs = element.GetStandardDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach (Tuple<XNode2D, StandardDOFType, int> entry in elementDofs)
            {
                int globalFreeDof;
                bool isFree = this.freeDofs.TryGetValue(entry.Item1, entry.Item2, out globalFreeDof);
                if (isFree) elementVector[entry.Item3] = globalFreeVector[globalFreeDof];
                else
                {
                    int globalConstrainedDof = this.constrainedDofs[entry.Item1, entry.Item2];
                    elementVector[entry.Item3] = globalConstrainedVector[globalConstrainedDof];
                }
            }
            return elementVector;
        }

        private static void EnumerateFreeDofs(IDictionary<XNode2D, HashSet<StandardDOFType>> nodalDOFTypes, 
            ITable<XNode2D, StandardDOFType, double> constraints, 
            out int freeDofsCount, out Table<XNode2D, StandardDOFType, int> freeDofs)
        {
            freeDofs = new Table<XNode2D, StandardDOFType, int>();
            int counter = 0;
            foreach (var pair in nodalDOFTypes)
            {
                XNode2D node = pair.Key;
                foreach (StandardDOFType dofType in pair.Value)
                {
                    if (!constraints.Contains(node, dofType)) freeDofs[node, dofType] = counter++;
                }
            }
            freeDofsCount = counter;
        }

        private static void EnumerateConstrainedDofs(IDictionary<XNode2D, HashSet<StandardDOFType>> nodalDOFTypes,
            ITable<XNode2D, StandardDOFType, double> constraints,
            out int constrainedDofsCount, out Table<XNode2D, StandardDOFType, int> constrainedDofs)
        {
            constrainedDofs = new Table<XNode2D, StandardDOFType, int>();
            int counter = 0;
            foreach (Tuple<XNode2D, StandardDOFType, double> entry in constraints)
            {
                constrainedDofs[entry.Item1, entry.Item2] = counter++;
            }
            constrainedDofsCount = counter;
        }

        private static IDictionary<XNode2D, HashSet<StandardDOFType>> FindUniqueDOFTypes(IEnumerable<XContinuumElement2D> elements)
        {
            var totalDofs = new SortedDictionary<XNode2D, HashSet<StandardDOFType>>();
            foreach (XContinuumElement2D element in elements)
            {
                ITable<XNode2D, StandardDOFType, int> elementDofs = element.GetStandardDofs();
                foreach (XNode2D node in elementDofs.GetRows())
                {
                    HashSet<StandardDOFType> dofsOfThisNode;
                    bool alreadyExists = totalDofs.TryGetValue(node, out dofsOfThisNode);
                    if (!alreadyExists)
                    {
                        dofsOfThisNode = new HashSet<StandardDOFType>();
                        totalDofs.Add(node, dofsOfThisNode);
                    }
                    dofsOfThisNode.UnionWith(elementDofs.GetColumnsOfRow(node));
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
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
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
