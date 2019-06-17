using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    internal static class SkylineSlicing
    {
        internal static double[] GetColumn(double[] skyValues, int[] skyDiagOffsets, int colIndex)
        {
            int order = skyDiagOffsets.Length - 1;
            var columnVector = new double[order];

            // Upper triangle and diagonal entries of the column are stored explicitly and contiguously
            int diagOffset = skyDiagOffsets[colIndex];
            int colHeight = skyDiagOffsets[colIndex + 1] - diagOffset - 1; // excluding diagonal
            for (int k = 0; k <= colHeight; ++k) columnVector[colIndex - k] = skyValues[diagOffset + k];

            // Lower triangle entries of the column can be found in the row with the same index
            for (int j = colIndex + 1; j < order; ++j)
            {
                int otherDiagOffset = skyDiagOffsets[j];
                int otherColHeight = skyDiagOffsets[j + 1] - otherDiagOffset - 1; // excluding diagonal
                int entryHeight = j - colIndex; // excluding diagonal
                if (entryHeight <= otherColHeight) // inside stored non zero pattern
                {
                    columnVector[j] = skyValues[otherDiagOffset + entryHeight];
                }
            }

            return columnVector;
        }

        internal static double[] GetDiagonal(double[] skyValues, int[] skyDiagOffsets)
        {
            int order = skyDiagOffsets.Length - 1;
            double[] diag = new double[order];
            for (int j = 0; j < order; ++j) diag[j] = skyValues[skyDiagOffsets[j]];
            return diag;
        }

        internal static CscMatrix GetSubmatrixCsc(double[] skyValues, int[] skyDiagOffsets, int[] rowsToKeep, int[] colsToKeep)
        {
            if ((rowsToKeep.Length == 0) || (colsToKeep.Length == 0))
            {
                return CscMatrix.CreateFromArrays(rowsToKeep.Length, colsToKeep.Length, new double[0], new int[0],
                    new int[1] { 0 }, false);
            }
            var submatrix = DokColMajor.CreateEmpty(rowsToKeep.Length, colsToKeep.Length);
            for (int subCol = 0; subCol < colsToKeep.Length; ++subCol)
            {
                int col = colsToKeep[subCol];
                int diagOffset = skyDiagOffsets[col];
                int colHeight = skyDiagOffsets[col + 1] - diagOffset - 1; // excluding diagonal
                for (int subRow = 0; subRow < rowsToKeep.Length; ++subRow)
                {
                    int row = rowsToKeep[subRow];
                    if (row <= col)
                    {
                        int entryHeight = col - row; // excluding diagonal
                        if (entryHeight <= colHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[diagOffset + entryHeight];
                            if (val != 0.0) submatrix[subRow, subCol] = val; // Skyline format stores many unnecessary zeros.
                        }
                    }
                    else // Transpose row <-> col. The cached column height and offset must be recalculated.
                    {
                        int transposedDiagOffset = skyDiagOffsets[row];
                        int transposedColHeight = skyDiagOffsets[row + 1] - transposedDiagOffset - 1; // excluding diagonal
                        int entryHeight = row - col; // excluding diagonal
                        if (entryHeight <= transposedColHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[transposedDiagOffset + entryHeight];
                            if (val != 0.0) submatrix[subRow, subCol] = val; // Skyline format stores many unnecessary zeros.
                        }
                    }
                }
            }
            return submatrix.BuildCscMatrix(true);
        }

        internal static Matrix GetSubmatrixSymmetricFull(double[] skyValues, int[] skyDiagOffsets, int[] rowsColsToKeep)
        {
            int subOrder = rowsColsToKeep.Length;
            var submatrix = Matrix.CreateZero(subOrder, subOrder);
            for (int subCol = 0; subCol < subOrder; ++subCol)
            {
                int col = rowsColsToKeep[subCol];
                int diagOffset = skyDiagOffsets[col];
                int colHeight = skyDiagOffsets[col + 1] - diagOffset - 1; // excluding diagonal
                for (int subRow = 0; subRow <= subCol; ++subRow)
                {
                    int row = rowsColsToKeep[subRow];
                    if (row <= col)
                    {
                        int entryHeight = col - row; // excluding diagonal
                        if (entryHeight <= colHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[diagOffset + entryHeight];
                            if (val != 0.0) // TODO: Not sure if this speeds up things
                            {
                                submatrix[subRow, subCol] = val;
                                submatrix[subCol, subRow] = val;
                            }
                        }
                    }
                    else // Transpose row <-> col. The cached column height and offset must be recalculated.
                    {
                        int transposedDiagOffset = skyDiagOffsets[row];
                        int transposedColHeight = skyDiagOffsets[row + 1] - transposedDiagOffset - 1; // excluding diagonal
                        int entryHeight = row - col; // excluding diagonal
                        if (entryHeight <= transposedColHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[transposedDiagOffset + entryHeight];
                            if (val != 0.0) // TODO: Not sure if this speeds up things
                            {
                                submatrix[subRow, subCol] = val;
                                submatrix[subCol, subRow] = val;
                            }
                        }
                    }
                }
            }
            return submatrix;
        }

        internal static SymmetricMatrix GetSubmatrixSymmetricPacked(double[] skyValues, int[] skyDiagOffsets, 
            int[] rowsColsToKeep)
        {
            int subOrder = rowsColsToKeep.Length;
            var submatrix = SymmetricMatrix.CreateZero(subOrder);
            for (int subCol = 0; subCol < subOrder; ++subCol)
            {
                int col = rowsColsToKeep[subCol];
                int diagOffset = skyDiagOffsets[col];
                int colHeight = skyDiagOffsets[col + 1] - diagOffset - 1; // excluding diagonal
                for (int subRow = 0; subRow <= subCol; ++subRow)
                {
                    int row = rowsColsToKeep[subRow];
                    if (row <= col)
                    {
                        int entryHeight = col - row; // excluding diagonal
                        if (entryHeight <= colHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[diagOffset + entryHeight];
                            if (val != 0.0) submatrix[subRow, subCol] = val;// TODO: Not sure if this speeds up things
                        }
                    }
                    else // Transpose row <-> col. The cached column height and offset must be recalculated.
                    {
                        int transposedDiagOffset = skyDiagOffsets[row];
                        int transposedColHeight = skyDiagOffsets[row + 1] - transposedDiagOffset - 1; // excluding diagonal
                        int entryHeight = row - col; // excluding diagonal
                        if (entryHeight <= transposedColHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[transposedDiagOffset + entryHeight];
                            if (val != 0.0) submatrix[subRow, subCol] = val;// TODO: Not sure if this speeds up things
                        }
                    }
                }
            }
            return submatrix;
        }

        //TODO: Move this method to SparsityPatternSymmetric, so that I can optimize access to its private data
        internal static SparsityPatternSymmetric GetSubmatrixSymmetricPattern(double[] skyValues, int[] skyDiagOffsets,
            int[] rowsColsToKeep)
        {
            //TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
            //      Schur complements more efficiently.

            int subOrder = rowsColsToKeep.Length;
            if (subOrder == 0) throw new NonMatchingDimensionsException("The original matrix had 0 rows and columns.");
            var submatrix = SparsityPatternSymmetric.CreateEmpty(subOrder);
            for (int subCol = 0; subCol < subOrder; ++subCol)
            {
                int col = rowsColsToKeep[subCol];
                int diagOffset = skyDiagOffsets[col];
                int colHeight = skyDiagOffsets[col + 1] - diagOffset - 1; // excluding diagonal
                for (int subRow = 0; subRow <= subCol; ++subRow)
                {
                    int row = rowsColsToKeep[subRow];
                    if (row <= col)
                    {
                        int entryHeight = col - row; // excluding diagonal
                        if (entryHeight <= colHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[diagOffset + entryHeight];
                            if (val != 0.0) submatrix.AddEntryUpper(subRow, subCol); // Skyline format stores many unnecessary zeros.
                        }
                    }
                    else // Transpose row <-> col. The cached column height and offset must be recalculated.
                    {
                        int transposedDiagOffset = skyDiagOffsets[row];
                        int transposedColHeight = skyDiagOffsets[row + 1] - transposedDiagOffset - 1; // excluding diagonal
                        int entryHeight = row - col; // excluding diagonal
                        if (entryHeight <= transposedColHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[transposedDiagOffset + entryHeight];
                            if (val != 0.0) submatrix.AddEntryUpper(subRow, subCol); // Skyline format stores many unnecessary zeros.
                        }
                    }
                }
            }
            return submatrix;
        }

        //TODO: Move this method to DokSymmetric, so that I can optimize access to its private data
        internal static SkylineMatrix GetSubmatrixSymmetricSkyline(double[] skyValues, int[] skyDiagOffsets, 
            int[] rowsColsToKeep) 
        {
            //TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
            //      Schur complements more efficiently.

            int subOrder = rowsColsToKeep.Length;
            if (subOrder == 0) return SkylineMatrix.CreateFromArrays(subOrder, new double[0], new int[1] { 0 }, false, false);
            var submatrix = DokSymmetric.CreateEmpty(subOrder);
            for (int subCol = 0; subCol < subOrder; ++subCol)
            {
                int col = rowsColsToKeep[subCol];
                int diagOffset = skyDiagOffsets[col];
                int colHeight = skyDiagOffsets[col + 1] - diagOffset - 1; // excluding diagonal
                for (int subRow = 0; subRow <= subCol; ++subRow)
                {
                    int row = rowsColsToKeep[subRow];
                    if (row <= col)
                    {
                        int entryHeight = col - row; // excluding diagonal
                        if (entryHeight <= colHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[diagOffset + entryHeight];
                            if (val != 0.0) submatrix.SetEntryUpper(subRow, subCol, val); // Skyline format stores many unnecessary zeros.
                        }
                    }
                    else // Transpose row <-> col. The cached column height and offset must be recalculated.
                    {
                        int transposedDiagOffset = skyDiagOffsets[row];
                        int transposedColHeight = skyDiagOffsets[row + 1] - transposedDiagOffset - 1; // excluding diagonal
                        int entryHeight = row - col; // excluding diagonal
                        if (entryHeight <= transposedColHeight) // inside stored non zero pattern
                        {
                            double val = skyValues[transposedDiagOffset + entryHeight];
                            if (val != 0.0) submatrix[subRow, subCol] = val; // Skyline format stores many unnecessary zeros.
                        }
                    }
                }
            }
            return submatrix.BuildSkylineMatrix();
        }

        internal static SkylineMatrix GetSubmatrixSymmetricSkylineLegacy(double[] skyValues, int[] skyDiagOffsets, 
            int[] rowsColsToKeep)
        {
            // Taken from G. Stavroulakis code in a FETI-DP solver.
            //TODO: It does not seem to work correctly if rowsColsToKeep is not sorted in ascending order.
            //TODO: compare against the version that uses DOK, in terms of a) time/memory required for creating the submatrix,
            //      b) quality of the submatrix the DOK version may avoid some zeros.

            int subOrder = rowsColsToKeep.Length;
            int[] subDiagOffsets = new int[subOrder + 1];
            int offset = 0;
            for (int i = 0; i < subOrder; i++)
            {
                subDiagOffsets[i] = offset;
                int col = rowsColsToKeep[i];
                int fromRow = col - skyDiagOffsets[col + 1] + skyDiagOffsets[col] + 1; // top non zero entry of this col
                var newRows = rowsColsToKeep.Count(x => x >= fromRow && x <= col);
                offset += newRows;
            }
            subDiagOffsets[subOrder] = offset;

            var subValues = new double[subDiagOffsets[subDiagOffsets.Length - 1]];
            offset = 0;
            for (int i = 0; i < subOrder; i++)
            {
                int col = rowsColsToKeep[i];
                int fromRow = col - skyDiagOffsets[col + 1] + skyDiagOffsets[col] + 1; // top non zero entry of this col
                var newRows = rowsColsToKeep.Where(x => x >= fromRow && x <= col).OrderByDescending(x => x).ToArray<int>();
                for (int j = 0; j < newRows.Length; j++)
                {
                    subValues[offset] = skyValues[skyDiagOffsets[col] + col - newRows[j]];
                    offset++;
                }
            }

            return SkylineMatrix.CreateFromArrays(subOrder, subValues, subDiagOffsets, false, false);
        }
    }
}
