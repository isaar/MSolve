using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// Mainly useful for checking other formats, since it is more difficult to mess up.
    /// </summary>
    static class GlobalDenseAssembler
    {
        public static (Matrix Kff, Matrix Kfc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            var Kff = Matrix.CreateZero(numDofsFree, numDofsFree);
            var Kfc = Matrix.CreateZero(numDofsFree, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrix and add its contributions to the global matrices.
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                AddElementToGlobalMatrix(Kff, kss, mapStandard, mapStandard);
                AddElementToGlobalMatrix(Kfc, kss, mapStandard, mapConstrained);

                // Build enriched element matrices and add their contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    (Matrix kee, Matrix kes) = element.BuildEnrichedStiffnessMatricesLower();
                    Matrix kse = kes.Transpose();
                    AddElementToGlobalMatrix(Kff, kes, mapEnriched, mapStandard);
                    AddElementToGlobalMatrix(Kff, kse, mapStandard, mapEnriched);
                    AddElementToGlobalMatrix(Kff, kee, mapEnriched, mapEnriched);
                    AddElementToGlobalMatrix(Kfc, kes, mapEnriched, mapConstrained);
                }
            }

            return (Kff, Kfc);
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
