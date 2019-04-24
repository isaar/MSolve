using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;

/// TODO: abstract this with an IGlobalEnumerator2D interface. That way I could have different dof orders:
/// e.g. 1) first all standard dofs and then all enriched, 2) First all dofs (std & enr) of the 1st node, 
/// then of the 2nd, etc.
/// TODO: The enumeration should be in a dedicated class DofOrderer. Managing and exposing the dof tables and 
/// the conversion methods should be in a DofManager class or in Model. The DofManager will greatly violate LoD 
/// (without wrapping access to the tables and most methods), while if they are in Model it will violate SRP
/// TODO: The MatchElementToGlobal... methods do a lot of duplicate iterations over the dofs. Perhaps a completely different way 
/// should be found, which also changes the way the GlobalAssembler works. The ideal would be to split each element 
/// matrix into the free, constrained and artificial parts required by the Assembler, without any overlaps, which 
/// are currently denoted as -1 and need to be checked every single time they are created and accessed.
namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    /// <summary>
    /// The enriched dofs are numbered after all standard
    /// </summary>
    class SeparateDofOrderer: DofOrdererBase
    {
        private SeparateDofOrderer(int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs,
            int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs,
            int standardDofsCount, DofTable<StructuralDof> standardDofs):
            base(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, standardDofsCount, standardDofs)
        {
        }

        public static SeparateDofOrderer Create(Model2D_old model)
        {
            // TODO: I should probably have a Constraint or Constraints class, to decouple this class from the collections 
            // used to represent constraints
            IDictionary<XNode, HashSet<StructuralDof>> nodalDofTypes = FindUniqueDofTypes(model.Elements);
            (int standardDofsCount, DofTable<StructuralDof> standardDofs) = 
                OrderStandardDofs(nodalDofTypes, model.Constraints);
            (int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs) =
                OrderConstrainedDofs(model.Constraints);
            (int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs) = 
                OrderEnrichedDofs(model.Nodes, standardDofsCount);

            #region DEBUG code
            //Console.WriteLine("------------------------ DEBUG ------------------------------");
            //Console.WriteLine("Free standard dofs: ");
            //Console.WriteLine(standardDofs);
            //Console.WriteLine("Enriched dofs: ");
            //Console.WriteLine(enrichedDofs);
            //Console.WriteLine("Constrained dofs: ");
            //Console.WriteLine(constrainedDofs);
            //Console.WriteLine("------------------------ /DEBUG ------------------------------");
            #endregion

            return new SeparateDofOrderer(constrainedDofsCount, constrainedDofs, enrichedDofsCount, enrichedDofs, 
                standardDofsCount, standardDofs);
        }

        // This is for problems that have rotational dofs only at some nodes. In continuum mechanics, we can just assign 2  
        // standard dofs per node and avoid iterating the elements. TODO: There should be different enumerators
        private static IDictionary<XNode, HashSet<StructuralDof>> FindUniqueDofTypes(IEnumerable<XContinuumElement2D> elements)
        {
            var totalDofs = new SortedDictionary<XNode, HashSet<StructuralDof>>();
            foreach (XContinuumElement2D element in elements)
            {
                ITable<XNode, StructuralDof, int> elementDofs = element.GetStandardDofs();
                foreach (XNode node in elementDofs.GetRows())
                {
                    bool alreadyExists = totalDofs.TryGetValue(node, out HashSet<StructuralDof> dofsOfThisNode);
                    if (!alreadyExists)
                    {
                        dofsOfThisNode = new HashSet<StructuralDof>();
                        totalDofs.Add(node, dofsOfThisNode);
                    }
                    dofsOfThisNode.UnionWith(elementDofs.GetColumnsOfRow(node));
                }
            }
            return totalDofs;
        }

        private static (int constrainedDofsCount, DofTable<StructuralDof> constrainedDofs) OrderConstrainedDofs(
            ITable<XNode, StructuralDof, double> constraints)
        {
            var constrainedDofs = new DofTable<StructuralDof>();
            int counter = 0;
            foreach ((XNode node, StructuralDof dofType, double displacement) in constraints)
            {
                constrainedDofs[node, dofType] = counter++;
            }
            return (counter, constrainedDofs);
        }

        // Each artificial dof has index that is node major, then enrichment item major, then enrichment function major and finally axis minor
        private static (int enrichedDofsCount, DofTable<EnrichedDof> enrichedDofs) OrderEnrichedDofs(
            IEnumerable<XNode> nodes, int standardDofsCount)
        {
            var enrichedDofs = new DofTable<EnrichedDof>();
            int dofCounter = standardDofsCount; // This if I put everything in the same matrix
            //int dofCounter = 0; // This if I use different matrices
            foreach (XNode node in nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    foreach (EnrichedDof dofType in enrichment.Dofs) // Are dofs determined by the element type (e.g. structural) as well?
                    {
                        enrichedDofs[node, dofType] = dofCounter++;
                    }
                }
            }
            return (dofCounter - standardDofsCount, enrichedDofs);
        }

        //Node major ordering
        private static (int standardDofsCount, DofTable<StructuralDof> standardDofs) OrderStandardDofs(
            IDictionary<XNode, HashSet<StructuralDof>> nodalDofTypes, ITable<XNode, StructuralDof, double> constraints)
        {
            var standardDofs = new DofTable<StructuralDof>();
            int counter = 0;
            foreach (var pair in nodalDofTypes)
            {
                XNode node = pair.Key;
                foreach (StructuralDof dofType in pair.Value)
                {
                    if (!constraints.Contains(node, dofType)) standardDofs[node, dofType] = counter++;
                }
            }
            return (counter, standardDofs);
        }
    }
}
