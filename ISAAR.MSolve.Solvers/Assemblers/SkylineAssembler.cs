using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

//TODO: The F = Ff - Kfc*Fc should not be done in the solver. The solver should only operate on the final linear systems.
//      It could be done here or in the analyzer.
//TODO: Optimizations should be available when building matrices with the same sparsity pattern. These optimizations should work
//      for all assemblers, not only skyline.
//TODO: checking whether the indexer is null is not enough. The assembler should rebuild the indexer whenever the dof ordering
//      changes.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is symmetric and stored in Skyline 
    /// format, which is suitable for the Cholesky factorization (e.g. in a direct solver).
    /// </summary>
    public class SkylineAssembler : IGlobalMatrixAssembler<SkylineMatrix>
    {
        private const string name = "SkylineAssembler"; // for error messages

        bool areIndexersCached = false;
        private SkylineBuilder skylineBuilder;
        private ConstrainedMatricesAssembler constrainedAssembler = new ConstrainedMatricesAssembler();

        public SkylineMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements,
            IElementMatrixProvider_v2 matrixProvider)
        {
            if (!areIndexersCached)
            {
                skylineBuilder = SkylineBuilder.Create(dofOrdering.NumFreeDofs,
                    FindSkylineColumnHeights(elements, dofOrdering.NumFreeDofs, dofOrdering.FreeDofs));
                areIndexersCached = true;
            }
            else skylineBuilder.ClearValues();

            foreach (IElement_v2 element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                skylineBuilder.AddSubmatrixSymmetric(elementMatrix, elementDofIndices, subdomainDofIndices);
            }

            return skylineBuilder.BuildSkylineMatrix();
        }

        public (SkylineMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree, 
            IMatrixView matrixConstrConstr) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering, 
            IEnumerable<IElement_v2> elements, IElementMatrixProvider_v2 matrixProvider)
        {
            if (!areIndexersCached)
            {
                skylineBuilder = SkylineBuilder.Create(freeDofOrdering.NumFreeDofs,
                    FindSkylineColumnHeights(elements, freeDofOrdering.NumFreeDofs, freeDofOrdering.FreeDofs));
                areIndexersCached = true;
            }
            else skylineBuilder.ClearValues();

            //TODO: also reuse the indexers of the constrained matrices.
            constrainedAssembler.InitializeNewMatrices(freeDofOrdering.NumFreeDofs, constrainedDofOrdering.NumConstrainedDofs);

            // Process the stiffness of each element
            foreach (IElement_v2 element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofsFree, int[] subdomainDofsFree) = freeDofOrdering.MapFreeDofsElementToSubdomain(element);
                (int[] elementDofsConstrained, int[] subdomainDofsConstrained) = 
                    constrainedDofOrdering.MapConstrainedDofsElementToSubdomain(element);

                IMatrix elementMatrix = matrixProvider.Matrix(element);
                skylineBuilder.AddSubmatrixSymmetric(elementMatrix, elementDofsFree, subdomainDofsFree);
                constrainedAssembler.AddElementMatrix(elementMatrix, elementDofsFree, subdomainDofsFree,
                    elementDofsConstrained, subdomainDofsConstrained);
            }

            // Create the free and constrained matrices. 
            SkylineMatrix matrixFreeFree = skylineBuilder.BuildSkylineMatrix();
            (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) = constrainedAssembler.BuildMatrices();
            return (matrixFreeFree, matrixConstrFree.TransposeToCSC(false), matrixConstrFree, matrixConstrConstr);
        }

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
            skylineBuilder = null;
            areIndexersCached = false;
        }

        //TODO: If one element engages some dofs (of a node) and another engages other dofs, the ones not in the intersection 
        // are not dependent from the rest. This method assumes dependency for all dofs of the same node. This is a rare occasion 
        // though.
        private static int[] FindSkylineColumnHeights(IEnumerable<IElement_v2> elements,
            int numFreeDofs, DofTable freeDofs)
        {
            int[] colHeights = new int[numFreeDofs]; //only entries above the diagonal count towards the column height
            foreach (IElement_v2 element in elements)
            {
                //TODO: perhaps I could use dofOrdering.MapFreeDofsElementToSubdomain(element). This way they can be cached,
                //      which would speed up the code when building the values array. However, if there is not enough memory for 
                //      caching, performance may take a hit since building the mapping arrays does redundant stuff (probably?).
                //      In any case, benchmarking is needed.
                //TODO: perhaps the 2 outer loops could be done at once to avoid a lot of dof indexing. Could I update minDof
                //      and colHeights[] at once? At least I could store the dofIndices somewhere

                IList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0.0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (var node in elementNodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (var node in elementNodes)
                {
                    foreach (int dof in freeDofs.GetValuesOfRow(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return colHeights;
        }
    }
}
