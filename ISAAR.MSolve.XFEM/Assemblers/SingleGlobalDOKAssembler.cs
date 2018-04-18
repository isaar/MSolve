using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

/// <summary>
/// The matrix that will be "inverted" is a unified DOK matrix, where enriched dofs are numbered after all 
/// standard dofs. 
/// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
/// needed. Also the global constrained submatrix must be sparse.
/// </summary>
namespace ISAAR.MSolve.XFEM.Assemblers
{
    static class SingleGlobalDOKAssembler
    {
        public static (DOKSymmetricColMajor Kff, Matrix Kfc) BuildGlobalMatrix(Model2D model)
        {
            DOFEnumerator dofEnumerator = model.DofEnumerator;
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;
            int allFreeDofsCount = dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            var Kff = DOKSymmetricColMajor.CreateEmpty(allFreeDofsCount);

            // TODO: this should be in a sparse format. Only used for SpMV and perhaps transpose SpMV.
            // Row = standard free dofs + enriched dofs. Columns = standard constrained dofs. 
            Matrix Kfc = Matrix.CreateZero(dofEnumerator.FreeDofsCount + dofEnumerator.ArtificialDofsCount,
                dofEnumerator.ConstrainedDofsCount);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                Matrix kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                // 2) The Elements have two versions of BuildStiffness(). 
                // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                Matrix kse = kes.Transpose(); 

                // Element to global dofs mappings
                // TODO: perhaps that could be done during the assembly to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element, 
                    out IReadOnlyDictionary<int, int> mapFree, 
                    out IReadOnlyDictionary<int, int> mapConstrained);
                IReadOnlyDictionary<int, int> mapEnriched = dofEnumerator.MatchElementToGlobalArtificialDofsOf(element);

                // Add the element contributions to the global matrices
                AddElementToGlobalMatrix(Kff, kss, mapFree, mapFree);
                AddElementToGlobalMatrix(Kff, kse, mapFree,  mapEnriched);
                AddElementToGlobalMatrix(Kff, kee, mapEnriched, mapEnriched);

                AddElementToGlobalMatrix(Kfc, kss, mapFree, mapConstrained);
                AddElementToGlobalMatrix(Kfc, kes, mapEnriched, mapConstrained);
            }

            // Test the global matrices
            //SingleGlobalSkylineAssembler.BuildGlobalMatrix(model, out SkylineMatrix2D KuuExpected, out Matrix2D KucExpected);
            //Console.WriteLine("Check Kuu: ");
            //CheckMatrix(Kuu, KuuExpected);
            //Console.WriteLine();
            //Console.WriteLine("Check Kuc: ");
            //CheckMatrix(Kuc, KucExpected);
            return (Kff, Kfc);
        }

        private static void AddElementToGlobalMatrix(DOKSymmetricColMajor globalMatrix, Matrix elementMatrix,
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

                    // SymmetricDOKColMajor matrix stores the upper triangle
                    // Check if the entry is in the upper triangle. This check concerns global dofs only.
                    // No need to check if the entry is within the row height.
                    // TODO: This check is redundant for the enriched-standard submatrices, since they always 
                    // have row=enrichedDof > col=standardDof. Perhaps I could have an optimized method that 
                    // doesn't check that height>0. Frankly it would be too much code just for avoiding this 
                    // check. Perhaps after optimizing other things first. Alternatively overload the method, with a version
                    // that accepts Symmetric matrices (where there would be a TraverseUpper perhaps)
                    if (globalCol >= globalRow)
                    {
                        globalMatrix.AddToEntry(globalRow, globalCol, elementMatrix[elementRow, elementCol]);
                    }
                }
            }
        }

        /// <summary>
        /// This method is for the non symmetric part of the matrix (free - constrained).
        /// TODO: The global matrix should be sparse. Probably in CSR/DOK format.
        /// </summary>
        /// <param name="globalMatrix"></param>
        /// <param name="elementMatrix"></param>
        /// <param name="elementRowsToGlobalRows"></param>
        /// <param name="elementColsToGlobalCols"></param>
        private static void AddElementToGlobalMatrix(Matrix globalMatrix, Matrix elementMatrix,
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

                    globalMatrix[globalRow, globalCol] = elementMatrix[elementRow, elementCol];
                }
            }
        }

        private static void CheckMatrix(Matrix computed, Matrix expected)
        {
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.NumRows != expected.NumRows) || (computed.NumColumns != expected.NumColumns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.NumRows; ++i)
            {
                for (int j = 0; j < computed.NumColumns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                    }
                }
            }
        }
    }
}
