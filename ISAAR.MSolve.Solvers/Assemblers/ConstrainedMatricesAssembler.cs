using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

//TODO: Reuse indexers. In that case they DOKs should probably be handled by the parent assembler. Therefore there is no reason
//      for this class.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Creates the subdomain level matrices that correspond to constrained dofs Kcf, Kcc, where f = free dof, c = constrained 
    /// dof, A * X = B, A = [ Aff Acf^T; Acf Acc ], X = [ Xf; Xc ], B = [ Bf, Bc ] (Matlab notation).
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ConstrainedMatricesAssembler
    {
        private DokRowMajor dokConstrFree, dokConstrConstr;

        internal void AddElementMatrix(IMatrixView elementMatrix, int[] elementDofsFree, int[] subdomainDofsFree, 
            int[] elementDofsConstrained, int[] subdomainDofsConstrained)
        {
            dokConstrFree.AddSubmatrix(elementMatrix, elementDofsConstrained, subdomainDofsConstrained, 
                elementDofsFree, subdomainDofsFree);
            dokConstrConstr.AddSubmatrixSymmetric(elementMatrix, elementDofsConstrained, subdomainDofsConstrained);
        }

        internal (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) BuildMatrices()
        {
            // Also free the memory used by each builders, as soon as it is no longer used.
            CsrMatrix matrixConstrFree = dokConstrFree.BuildCsrMatrix(true);
            dokConstrFree = null;
            CsrMatrix matrixConstrConstr = dokConstrConstr.BuildCsrMatrix(true);
            dokConstrConstr = null;
            return (matrixConstrFree, matrixConstrConstr);
        }

        internal void InitializeNewMatrices(int numFreeDofs, int numConstrainedDofs)
        {
            dokConstrFree = DokRowMajor.CreateEmpty(numConstrainedDofs, numFreeDofs);
            dokConstrConstr = DokRowMajor.CreateEmpty(numConstrainedDofs, numConstrainedDofs);
        }
    }
}
