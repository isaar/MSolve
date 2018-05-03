using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    interface IDOFEnumerator
    {
        int ConstrainedDofsCount { get; }
        int EnrichedDofsCount { get ; }
        int FreeDofsCount { get; }

        IDOFEnumerator DeepCopy();

        int GetFreeDofOf(XNode2D node, DisplacementDOF dofType);
        IEnumerable<int> GetFreeDofsOf(XNode2D node);
        List<int> GetFreeDofsOf(XContinuumElement2D element);

        int GetConstrainedDofOf(XNode2D node, DisplacementDOF dofType);
        IEnumerable<int> GetConstrainedDofsOf(XNode2D node);
        List<int> GetConstrainedDofsOf(XContinuumElement2D element); // Also add a method that simultaneously returns free+constrained

        int GetEnrichedDofOf(XNode2D node, EnrichedDOF dofType);
        IEnumerable<int> GetEnrichedDofsOf(XNode2D node);
        List<int> GetEnrichedDofsOf(XContinuumElement2D element);

        void MatchElementToGlobalStandardDofsOf(XContinuumElement2D element,
            out IReadOnlyDictionary<int, int> elementToGlobalFreeDofs,
            out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs);

        IReadOnlyDictionary<int, int> MatchElementToGlobalEnrichedDofsOf(XContinuumElement2D element);

        Vector ExtractDisplacementVectorOfElementFromGlobal(XContinuumElement2D element,
            Vector globalFreeVector, Vector globalConstrainedVector);

        Vector ExtractEnrichedDisplacementsOfElementFromGlobal(XContinuumElement2D element, Vector globalFreeVector);

        double[,] GatherNodalDisplacements(Model2D model, Vector solution);

        ITable<XNode2D, EnrichedDOF, double> GatherEnrichedNodalDisplacements(Model2D model, Vector solution);

        /// <summary>
        /// Renumbers the dof indices according th the given permutation vector and direction. 
        /// If (<paramref name="oldToNew"/> == true), then newIndex[dof] = <paramref name="permutation"/>[oldIndex[dof]].
        /// Else oldIndex[dof] = <paramref name="permutation"/>[nwIndex[dof]]
        /// </summary>
        /// <param name="permutation">The permutation vector.</param>
        /// <param name="oldToNew">The direction it should be applied to.</param>
        void ReorderUnconstrainedDofs(IReadOnlyList<int> permutation, bool oldToNew);

        void WriteToConsole();
    }
}
