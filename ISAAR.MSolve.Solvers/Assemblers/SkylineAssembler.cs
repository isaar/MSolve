using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.Ordering;

//TODO: The F = Ff - Kfc*Fc should not be done in the solver. The solver should only operate on the final linear systems.
//      It could be done here or in the analyzer.
//TODO: Optimizations should be available when building matrices with the same sparsity pattern. These optimizations should work
//      for all assemblers, not only skyline.
//TODO: remove casts and the logic for specific elements
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineAssembler
    {
        public (SkylineMatrix Kff, DokRowMajor Kfc) BuildGlobalMatrices(IEnumerable<IElement> elements,
            AllDofOrderer dofOrderer, IElementMatrixProvider matrixProvider)
        {
            int numConstrainedDofs = dofOrderer.NumConstrainedDofs;
            int numFreeDofs = dofOrderer.NumFreeDofs;
            SkylineBuilder Kff = FindSkylineColumnHeights(elements, numFreeDofs, dofOrderer.FreeDofs);
            var Kfc = DokRowMajor.CreateEmpty(numFreeDofs, numConstrainedDofs);

            foreach (IElement elementWrapper in elements)
            {
                ContinuumElement2D element = (ContinuumElement2D)(elementWrapper.IElementType);
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (IReadOnlyDictionary<int, int> mapStandard, IReadOnlyDictionary<int, int> mapConstrained) =
                    dofOrderer.MapDofsElementToGlobal(element);
                Matrix k = Conversions.MatrixOldToNew(matrixProvider.Matrix(elementWrapper));
                Kff.AddSubmatrixSymmetric(k, mapStandard);
                Kfc.AddSubmatrix(k, mapStandard, mapConstrained);
            }

            //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. dok.Build())
            return (Kff.BuildSkylineMatrix(), Kfc);
        }

        public SkylineMatrix BuildGlobalMatrix(IEnumerable<IElement> elements, FreeDofOrderer dofOrderer, 
            IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = dofOrderer.NumFreeDofs;
            SkylineBuilder Kff = FindSkylineColumnHeights(elements, numFreeDofs, dofOrderer.FreeDofs);

            foreach (IElement elementWrapper in elements)
            {
                ContinuumElement2D element = (ContinuumElement2D)(elementWrapper.IElementType);
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                IReadOnlyDictionary<int, int> mapStandard  = dofOrderer.MapFreeDofsElementToGlobal(element);
                Matrix k = Conversions.MatrixOldToNew(matrixProvider.Matrix(elementWrapper));
                Kff.AddSubmatrixSymmetric(k, mapStandard);
            }

            return Kff.BuildSkylineMatrix();
        }

        //TODO: If one element engages some dofs (of a node) and another engages other dofs, the ones not in the intersection 
        // are not dependent from the rest. This method assumes dependency for all dofs of the same node. This is a rare occasion 
        // though.
        private static SkylineBuilder FindSkylineColumnHeights(IEnumerable<IElement> elements,
            int numFreeDofs, DofTable<IDof> freeDofs)
        {
            int[] colHeights = new int[numFreeDofs]; //only entries above the diagonal count towards the column height
            foreach (IElement elementWrapper in elements)
            {
                ContinuumElement2D element = (ContinuumElement2D)(elementWrapper.IElementType);
                //TODO: perhaps the 2 outer loops could be done at once to avoid a lot of dof indexing. Could I update minDof
                //      and colHeights[] at once?

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (var node in element.Nodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (var node in element.Nodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                    foreach (int dof in freeDofs.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return SkylineBuilder.Create(numFreeDofs, colHeights);
        }
    }
}
