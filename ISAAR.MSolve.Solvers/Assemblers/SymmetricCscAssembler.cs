using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is in symmetric CSC format, namely only 
    /// the upper triangle is explicitly stored.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricCscAssembler : IGlobalMatrixAssembler<SymmetricCscMatrix>
    {
        private const string name = "SymmetricCscAssembler"; // for error messages
        private readonly bool sortColsOfEachRow;

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

        public SymmetricCscMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var Kff = DokSymmetric.CreateEmpty(numFreeDofs);

            foreach (IElement element in elements)
            {
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                Matrix k = matrixProvider.Matrix(element).LegacyToNewMatrix();
                Kff.AddSubmatrixSymmetric(k, elementDofIndices, subdomainDofIndices);
            }

            return Kff.BuildSymmetricCscMatrix(sortColsOfEachRow);
        }

        public SymmetricCscMatrix BuildGlobalMatrix(ISubdomain subdomain, IElementMatrixProvider matrixProvider)
        {
            throw new NotImplementedException();
        }
    }
}
