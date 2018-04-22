using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Utilities;

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
namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    /// <summary>
    /// The enriched dofs are numbered after all standard
    /// </summary>
    class DOFEnumeratorSeparate: DOFEnumeratorBase
    {
        private DOFEnumeratorSeparate(int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs,
            int enrichedDofsCount, Table<XNode2D, EnrichedDOF, int> enrichedDofs,
            int freeDofsCount, Table<XNode2D, DisplacementDOF, int> freeDofs):
            base(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, freeDofsCount, freeDofs)
        {
        }

        public static DOFEnumeratorSeparate Create(Model2D model)
        {
            // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections 
            // used to represent constraints
            IDictionary<XNode2D, HashSet<DisplacementDOF>> nodalDOFTypes = FindUniqueDOFTypes(model.Elements);
            (int freeDofsCount, Table< XNode2D, DisplacementDOF, int> freeDofs) = 
                EnumerateFreeDofs(nodalDOFTypes, model.Constraints);
            (int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs) =
                EnumerateConstrainedDofs(model.Constraints);
            (int enrichedDofsCount, Table<XNode2D, EnrichedDOF, int> enrichedDofs) = 
                EnumerateEnrichedDofs(model.Nodes, freeDofsCount);
            return new DOFEnumeratorSeparate(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, 
                freeDofsCount, freeDofs);
        }

        private static (int constrainedDofsCount, Table<XNode2D, DisplacementDOF, int> constrainedDofs) EnumerateConstrainedDofs(
            ITable<XNode2D, DisplacementDOF, double> constraints)
        {
            var constrainedDofs = new Table<XNode2D, DisplacementDOF, int>();
            int counter = 0;
            foreach (Tuple<XNode2D, DisplacementDOF, double> entry in constraints)
            {
                constrainedDofs[entry.Item1, entry.Item2] = counter++;
            }
            return (counter, constrainedDofs);
        }

        // Each artificial dof has index that is node major, then enrichment item major, then enrichment function major and finally axis minor
        private static (int enrichedDofsCount, Table<XNode2D, EnrichedDOF, int> enrichedDofs) EnumerateEnrichedDofs(
            IEnumerable<XNode2D> nodes, int standardDofsCount)
        {
            var enrichedDofs = new Table<XNode2D, EnrichedDOF, int>();
            int dofCounter = standardDofsCount; // This if I put everything in the same matrix
            //int dofCounter = 0; // This if I use different matrices
            foreach (XNode2D node in nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDOF dofType in enrichment.DOFs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        enrichedDofs[node, dofType] = dofCounter++;
                    }
                }
            }
            return (dofCounter - standardDofsCount, enrichedDofs);
        }

        //Node major ordering
        private static (int freeDofsCount, Table<XNode2D, DisplacementDOF, int> freeDofs) EnumerateFreeDofs(
            IDictionary<XNode2D, HashSet<DisplacementDOF>> nodalDOFTypes, ITable<XNode2D, DisplacementDOF, double> constraints)
        {
            var freeDofs = new Table<XNode2D, DisplacementDOF, int>();
            int counter = 0;
            foreach (var pair in nodalDOFTypes)
            {
                XNode2D node = pair.Key;
                foreach (DisplacementDOF dofType in pair.Value)
                {
                    if (!constraints.Contains(node, dofType)) freeDofs[node, dofType] = counter++;
                }
            }
            return (counter, freeDofs);
        }

        // This is for problems that have rotational dofs only at some nodes. In continuum mechanics, we can just assign 2  
        // standard dofs per node and avoid iterating the elements.
        private static IDictionary<XNode2D, HashSet<DisplacementDOF>> FindUniqueDOFTypes(IEnumerable<XContinuumElement2D> elements)
        {
            var totalDofs = new SortedDictionary<XNode2D, HashSet<DisplacementDOF>>();
            foreach (XContinuumElement2D element in elements)
            {
                ITable<XNode2D, DisplacementDOF, int> elementDofs = element.GetStandardDofs();
                foreach (XNode2D node in elementDofs.GetRows())
                {
                    bool alreadyExists = totalDofs.TryGetValue(node, out HashSet<DisplacementDOF> dofsOfThisNode);
                    if (!alreadyExists)
                    {
                        dofsOfThisNode = new HashSet<DisplacementDOF>();
                        totalDofs.Add(node, dofsOfThisNode);
                    }
                    dofsOfThisNode.UnionWith(elementDofs.GetColumnsOfRow(node));
                }
            }
            return totalDofs;
        }
    }
}
