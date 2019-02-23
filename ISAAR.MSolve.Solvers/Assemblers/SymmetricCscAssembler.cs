using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;

//TODO: Instead of storing the raw CSC arrays, use a reusable DOK or SymmCscIndexer class. That class should provide methods to 
//      assemble the values part of the global matrix more efficiently than the general purpose DOK. The general purpose DOK 
//      should only be used to assemble the first global matrix and whenever the dof ordering changes. Now it is used everytime 
//      and the indexing arrays are discarded.
//TODO: I could also cache the symbolic factorization of SuiteSparse and reuse it. That would really speed up things.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in symmetric CSC format, namely only 
    /// the upper triangle is explicitly stored. This format is suitable for the SuiteSparse library and solvers that use it.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricCscAssembler : IGlobalMatrixAssembler<SymmetricCscMatrix>
    {
        private const string name = "SymmetricCscAssembler"; // for error messages
        private readonly bool sortColsOfEachRow;
        private ConstrainedMatricesAssembler constrainedAssembler = new ConstrainedMatricesAssembler();

        bool isIndexerCached = false;
        private int[] cachedRowIndices, cachedColOffsets;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sortColsOfEachRow">
        /// Sorting the columns of each row in the CSC storage format may increase performance of the factorization and 
        /// back/forward substitutions. It is recommended to set it to true.
        /// </param>
        public SymmetricCscAssembler(bool sortColsOfEachRow = true)
        {
            this.sortColsOfEachRow = sortColsOfEachRow;
        }

        public SymmetricCscMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement_v2> elements,
            IElementMatrixProvider_v2 matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var subdomainMatrix = DokSymmetric.CreateEmpty(numFreeDofs);

            foreach (IElement_v2 element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofIndices, subdomainDofIndices);
            }

            (double[] values, int[] rowIndices, int[] colOffsets) = subdomainMatrix.BuildSymmetricCscArrays(sortColsOfEachRow);
            if (!isIndexerCached)
            {
                cachedRowIndices = rowIndices;
                cachedColOffsets = colOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedRowIndices, rowIndices));
                Debug.Assert(Utilities.AreEqual(cachedColOffsets, colOffsets));
            }


            return SymmetricCscMatrix.CreateFromArrays(numFreeDofs, values, cachedRowIndices, cachedColOffsets, false);
        }

        public (SymmetricCscMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement_v2> elements, IElementMatrixProvider_v2 matrixProvider)
        {
            int numFreeDofs = freeDofOrdering.NumFreeDofs;
            var subdomainMatrix = DokSymmetric.CreateEmpty(numFreeDofs);

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
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofsFree, subdomainDofsFree);
                constrainedAssembler.AddElementMatrix(elementMatrix, elementDofsFree, subdomainDofsFree,
                    elementDofsConstrained, subdomainDofsConstrained);
            }

            // Create and cache the CSC arrays for the free dofs.
            (double[] values, int[] rowIndices, int[] colOffsets) = subdomainMatrix.BuildSymmetricCscArrays(sortColsOfEachRow);
            if (!isIndexerCached)
            {
                cachedRowIndices = rowIndices;
                cachedColOffsets = colOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedRowIndices, rowIndices));
                Debug.Assert(Utilities.AreEqual(cachedColOffsets, colOffsets));
            }

            // Create the free and constrained matrices. 
            subdomainMatrix = null; // Let the DOK be garbaged collected early, in case there isn't sufficient memory.
            var matrixFreeFree =
                SymmetricCscMatrix.CreateFromArrays(numFreeDofs, values, cachedRowIndices, cachedColOffsets, false);
            (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) = constrainedAssembler.BuildMatrices();
            return (matrixFreeFree, matrixConstrFree.TransposeToCSC(false), matrixConstrFree, matrixConstrConstr);
        }

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
            cachedRowIndices = null;
            cachedColOffsets = null;
            isIndexerCached = false;
        }
    }
}
