using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

//TODO: It has bugs. Stop using this thing...
namespace ISAAR.MSolve.XFEM.Assemblers
{
    /// <summary>
    /// The matrix that will be "inverted" is a unified Skyline matrix, where enriched dofs are numbered after all 
    /// standard dofs. 
    /// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
    /// needed. Also the global constrained submatrix must be sparse.
    /// </summary>
    static class SingleGlobalSkylineAssemblerOLD
    {
        public static (SkylineMatrix2D Kff, Matrix Kfc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            SkylineMatrix2D Kuu = new SkylineMatrix2D(CalculateSkylineIndexer(model, dofEnumerator));

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
                    AddElementToGlobalMatrix(Kuu, kes, elementToGlobalEnrDofs, elementToGlobalFreeDofs);
                    AddElementToGlobalMatrix(Kuu, kee, elementToGlobalEnrDofs, elementToGlobalEnrDofs);
                    AddElementToGlobalMatrix(Kuc, kes, elementToGlobalEnrDofs, elementToGlobalConstrainedDofs);
                }
            }

            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofEnumerator);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            return (Kuu, Kuc);
        }

        private static int[] CalculateSkylineIndexer(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int[] rowHeights = CalculateRowHeights(model, dofEnumerator);

            int[] rowIndex = new int[rowHeights.Length + 1]; // The indexer uses 1 extra ghost entry to facilitate some operations
            rowIndex[0] = 0; // The diagonal entries of the first 2 rows always have the same positions
            rowIndex[1] = 1;
            for (int i = 1; i < rowHeights.Length; ++i)
                rowIndex[i + 1] = rowIndex[i] + rowHeights[i] + 1; //+1 because the heights don't include the diagonal
            return rowIndex;
        }

        /// Only the entries above the diagonal are considered
        private static int[] CalculateRowHeights(Model2D model, IDOFEnumerator dofEnumerator) // I vote to call them RowWidths!!!
        {
            int[] rowHeights = new int[dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount];
            foreach (XContinuumElement2D element in model.Elements)
            {
                // To determine the row height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDOF = Int32.MaxValue;
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node)) minDOF = Math.Min(dof, minDOF);
                    foreach (int dof in dofEnumerator.GetEnrichedDofsOf(node)) Math.Min(dof, minDOF);
                }

                // The height of each row is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node))
                    {
                        rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                    foreach (int dof in dofEnumerator.GetEnrichedDofsOf(node))
                    {
                        rowHeights[dof] = Math.Max(rowHeights[dof], dof - minDOF);
                    }
                }
            }
            return rowHeights;
        }

        private static void AddElementToGlobalMatrix(SkylineMatrix2D globalMatrix, IMatrixView elementMatrix,
            IReadOnlyDictionary<int, int> elementRowsToGlobalRows, IReadOnlyDictionary<int, int> elementColsToGlobalCols)
        {
            foreach (var rowPair in elementRowsToGlobalRows)
            {
                int elementRow = rowPair.Key;
                foreach (var colPair in elementColsToGlobalCols)
                {
                    int elementCol = colPair.Key;
                    int globalRow, globalCol;
                    if (rowPair.Value >= colPair.Value)
                    {
                        globalRow = rowPair.Value;
                        globalCol = colPair.Value;
                    }
                    else
                    {
                        globalRow = colPair.Value;
                        globalCol = rowPair.Value;
                    }

                    int height = globalRow - globalCol; // Skyline matrix stores the lower triangle
                    // Check if the entry is in the lower triangle. This check concerns global dofs only.
                    // No need to check if the entry is within the row height.
                    // TODO: This check is redundant for the enriched-standard submatrices, since they always 
                    // have row=enrichedDof > col=standardDof. Perhaps I could have an optimized method that 
                    // doesn't check that height>0. Frankly it would be too much code just for avoiding this 
                    // check. Perhaps after optimizing other things first.
                    if (height >= 0) //this is probably not needed
                    {
                        // Shouldn't the skyline indexing logic be abstracted?
                        int skylinePos = globalMatrix.RowIndex[globalRow] + height;
                        globalMatrix.Data[skylinePos] += elementMatrix[elementRow, elementCol];
                    }
                }
            }
        }

        private static void AddSymmetricElementToGlobalMatrix(SkylineMatrix2D globalMatrix, IMatrixView elementMatrix,
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
                    int height = globalRow - globalCol; // Skyline matrix stores the lower triangle
                    // Check if the entry is in the lower triangle. This check concerns global dofs only.
                    // No need to check if the entry is within the row height.
                    // TODO: This check is redundant for the enriched-standard submatrices, since they always 
                    // have row=enrichedDof > col=standardDof. Perhaps I could have an optimized method that 
                    // doesn't check that height>0. Frankly it would be too much code just for avoiding this 
                    // check. Perhaps after optimizing other things first.
                    if (height >= 0)
                    {
                        // Shouldn't the skyline indexing logic be abstracted?
                        int skylinePos = globalMatrix.RowIndex[globalRow] + height;
                        globalMatrix.Data[skylinePos] += elementMatrix[elementRow, elementCol];
                    }
                }
            }
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

        private static void CheckMatrix(IIndexable2D expected, IIndexable2D computed)
        {
            bool isCorrect = true;
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
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Correct");
        }

        private static void CheckMatrix(IIndexable2D expected, IMatrix2D computed)
        {
            bool isCorrect = true;
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.Rows != expected.NumRows) || (computed.Columns != expected.NumColumns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.Rows; ++i)
            {
                for (int j = 0; j < computed.Columns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Correct");
        }
    }
}
