using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

//TODO: Clean this up: Not all these methods are actually necessary. Some are only used for testing specific orderers.  
namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    interface IDofOrderer
    {
        int NumConstrainedDofs { get; }
        int NumEnrichedDofs { get ; }
        int NumStandardDofs { get; }

        IDofOrderer DeepCopy();

        Vector ExtractDisplacementVectorOfElementFromGlobal(XContinuumElement2D element,
            Vector globalFreeVector, Vector globalConstrainedVector);

        Vector ExtractEnrichedDisplacementsOfElementFromGlobal(XContinuumElement2D element, Vector globalFreeVector);

        double[,] GatherNodalDisplacements(Model2D model, Vector solution);

        ITable<XNode2D, EnrichedDof, double> GatherEnrichedNodalDisplacements(Model2D model, Vector solution);

        int GetConstrainedDofOf(XNode2D node, DisplacementDof dofType);
        IEnumerable<int> GetConstrainedDofsOf(XNode2D node);
        List<int> GetConstrainedDofsOf(XContinuumElement2D element); // Also add a method that simultaneously returns free+constrained

        int GetEnrichedDofOf(XNode2D node, EnrichedDof dofType);
        IEnumerable<int> GetEnrichedDofsOf(XNode2D node);
        List<int> GetEnrichedDofsOf(XContinuumElement2D element);

        int GetStandardDofOf(XNode2D node, DisplacementDof dofType);
        IEnumerable<int> GetStandardDofsOf(XNode2D node);
        List<int> GetStandardDofsOf(XContinuumElement2D element);

        void MatchElementToGlobalStandardDofsOf(XContinuumElement2D element,
            out IReadOnlyDictionary<int, int> elementToGlobalStandardDofs,
            out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs);

        IReadOnlyDictionary<int, int> MatchElementToGlobalEnrichedDofsOf(XContinuumElement2D element);

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
