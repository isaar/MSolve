using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// Mainly useful for checking other formats, since it is more difficult to mess up.
    /// </summary>
    static class GlobalDenseAssembler
    {
        public static (Matrix Kff, Matrix Kfc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;
            int freeEnrichedDofsCount = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            Matrix Kuu = Matrix.CreateZero(freeEnrichedDofsCount, freeEnrichedDofsCount);

            // TODO: this should be in a sparse format. Only used for SpMV and perhaps transpose SpMV.
            // Row = standard free dofs + enriched dofs. Columns = standard constrained dofs. 
            Matrix Kuc = Matrix.CreateZero(dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount,
                dofEnumerator.ConstrainedDofsCount);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrices and add it contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> elementToGlobalFreeDofs,
                    out IReadOnlyDictionary<int, int> elementToGlobalConstrainedDofs);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                AddElementToGlobalMatrix(Kuu, kss, elementToGlobalFreeDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(Kuc, kss, elementToGlobalFreeDofs, elementToGlobalConstrainedDofs);

                // Build enriched element matrices and add it contributions to the global matrices
                IReadOnlyDictionary<int, int> elementToGlobalEnrDofs =
                    dofEnumerator.MatchElementToGlobalEnrichedDofsOf(element);
                if (elementToGlobalEnrDofs.Count > 0)
                {
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                    // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                    // 2) The Elements have two versions of BuildStiffness(). 
                    // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                    Matrix kse = kes.Transpose();
                    AddElementToGlobalMatrix(Kuu, kes, elementToGlobalEnrDofs, elementToGlobalFreeDofs);
                    AddElementToGlobalMatrix(Kuu, kes.Transpose(), elementToGlobalFreeDofs, elementToGlobalEnrDofs);
                    AddElementToGlobalMatrix(Kuu, kee, elementToGlobalEnrDofs, elementToGlobalEnrDofs);
                    AddElementToGlobalMatrix(Kuc, kes, elementToGlobalEnrDofs, elementToGlobalConstrainedDofs);
                }
            }

            return (Kuu, Kuc);
        }

        // TODO: The global matrix should be sparse. Probably in CSR/DOK format
        private static void AddElementToGlobalMatrix(Matrix globalMatrix, IMatrixView elementMatrix,
            IReadOnlyDictionary<int, int> elementRowsToGlobalRows, IReadOnlyDictionary<int, int> elementColsToGlobalCols)
        {
            foreach (var rowPair in elementRowsToGlobalRows)
            {
                int elementRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in elementColsToGlobalCols)
                {
                    int elementCol = colPair.Key;
                    int globalCol = colPair.Value;

                    globalMatrix[globalRow, globalCol] += elementMatrix[elementRow, elementCol];
                }
            }
        }
    }
}
