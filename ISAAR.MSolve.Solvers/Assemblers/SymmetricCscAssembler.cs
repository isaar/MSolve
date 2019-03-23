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
                IMatrix elementK = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementK, elementDofIndices, subdomainDofIndices);
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

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
            cachedRowIndices = null;
            cachedColOffsets = null;
            isIndexerCached = false;
        }
    }
}
