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
    static class DenseGlobalAssembler
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
                // Element matrices
                Matrix kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);


                // Element to global dofs mappings
                // TODO: perhaps that could be done during the assembly to avoid iterating over the dofs twice
                IReadOnlyDictionary<int, int> elementToGlobalFreeDofs, elementToGlobalConstrainedDofs;
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element,
                    out elementToGlobalFreeDofs, out elementToGlobalConstrainedDofs);
                IReadOnlyDictionary<int, int> elementToGlobalEnrDofs =
                    dofEnumerator.MatchElementToGlobalEnrichedDofsOf(element);

                // Add the element contributions to the global matrices
                AddElementToGlobalMatrix(Kuu, kss, elementToGlobalFreeDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(Kuu, kes, elementToGlobalEnrDofs, elementToGlobalFreeDofs);
                AddElementToGlobalMatrix(Kuu, kes.Transpose(), elementToGlobalFreeDofs, elementToGlobalEnrDofs);
                AddElementToGlobalMatrix(Kuu, kee, elementToGlobalEnrDofs, elementToGlobalEnrDofs);

                AddElementToGlobalMatrix(Kuc, kss, elementToGlobalFreeDofs, elementToGlobalConstrainedDofs);
                AddElementToGlobalMatrix(Kuc, kes, elementToGlobalEnrDofs, elementToGlobalConstrainedDofs);
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
