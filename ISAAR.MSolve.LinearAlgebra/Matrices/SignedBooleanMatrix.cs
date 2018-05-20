using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Each row corresponds to a displacement continuity equation. If this matrix is for the whole domain, sll rows will have 
    /// exactly 2 non zero entries (1 and -1). If it is only for a subdomain then some rows may be empty.
    /// Each column corresponds to a dof of one of the subdomains. Columns that correspond to dofs with multiplicity = 2 will 
    /// have 1 non zero entry (1 or -1). Columns that correspond to dofs with multiplicity > 2, will have at most 2 non zero 
    /// entries (1 and/or -1). The other columns do not correspond to boundary dofs and will be empty.
    /// </summary>
    public class SignedBooleanMatrix: IIndexable2D, ISparseMatrix
    {
        /// <summary>
        /// (row, (column, sign))
        /// </summary>
        private readonly Dictionary<int, Dictionary<int, int>> data;

        public SignedBooleanMatrix(int numRows, int numColumns)
        {
            this.NumRows = numRows;
            this.NumColumns = numColumns;
            this.data = new Dictionary<int, Dictionary<int, int>>();
        }

        public int NumRows { get; }
        public int NumColumns { get; }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                if (data.TryGetValue(rowIdx, out Dictionary<int, int> colSigns))
                {
                    if (colSigns.TryGetValue(colIdx, out int sign)) return sign;
                }
                return 0.0;
            }
        }

        /// <summary>
        /// Sets B[<paramref name="rowIdx"/>, <paramref name="colIdx"/>] = 1 or -1.
        /// </summary>
        /// <param name="rowIdx">The row index. This row must be empty before a non zero entry is added.</param>
        /// <param name="colIdx">The column index for the row <paramref name="rowIdx"/>.</param>
        /// <param name="sign">True for +1, false for -1</param>
        public void AddEntry(int rowIdx, int colIdx, bool sign)
        {
            if (data.TryGetValue(rowIdx, out Dictionary<int, int> colSigns))
            {
                colSigns.Add(colIdx, (sign ? 1 : -1));
            }
            else
            {
                var newColSigns = new Dictionary<int, int>();
                newColSigns.Add(colIdx, (sign ? 1 : -1));
                data.Add(rowIdx, newColSigns);
            }
        }

        //TODO: perhaps I should have a dedicated builder class. 
        //      Then the checks could also be more focused depending on if the matrix is global or for a subdomain
        public void CheckMatrix(bool forSubdomain) 
        {
            //TODO: dedicated exception class
            foreach (var wholeRow in data)
            {
                int[] signs = wholeRow.Value.Values.ToArray();
                if (signs.Length > 2) throw new MatrixPatternException(
                    $"Each row may have at most 2 entries, but row {wholeRow.Key} has {signs.Length}");
                else if ((signs.Length == 2) && (signs[0] != -signs[1])) throw new MatrixPatternException(
                    $"Row {wholeRow.Key} must have two opposite signs, but they were {signs[0]} and {signs[1]}");
            }
        }

        public double[,] CopyToArray2D()
        {
            return DenseStrategies.CopyToArray2D(this);
        }

        public Matrix CopyToFullMatrix(bool transpose)
        {
            // TODO: perhaps I should work with th col major arrays.
            if (transpose)
            {
                var dense = Matrix.CreateZero(this.NumColumns, this.NumRows);
                foreach (var wholeRow in data)
                {
                    foreach (var colValuePair in wholeRow.Value) dense[colValuePair.Key, wholeRow.Key] = colValuePair.Value;
                }
                return dense;
            }
            else
            {
                var dense = Matrix.CreateZero(this.NumRows, this.NumColumns);
                foreach (var wholeRow in data)
                {
                    foreach (var colValuePair in wholeRow.Value) dense[wholeRow.Key, colValuePair.Key] = colValuePair.Value;
                }
                return dense;
            }
        }

        public int CountNonZeros()
        {
            int count = 0;
            foreach (var wholeRow in data.Values) count += wholeRow.Count;
            return count;
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            foreach (var wholeRow in data)
            {
                foreach (var colVal in wholeRow.Value)
                {
                    yield return (wholeRow.Key, colVal.Key, colVal.Value);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public SparseFormat GetSparseFormat()
        {
            throw new NotImplementedException();
        }

        //TODO: I think that dealing with arrays will be faster than iterating the dictionaries. Another reason to separate 
        //      construction from multiplications.
        public Vector MultiplyRight(Vector vector, bool transposeThis)
        {
            if (transposeThis) return MultiplyTransposed(vector);
            else return MultiplyUntransposed(vector);
        }

        /// <summary>
        /// Not efficient. Meant for testing purposes.
        /// </summary>
        /// <returns></returns>
        public SignedBooleanMatrix Transpose()
        {
            var transpose = new SignedBooleanMatrix(NumColumns, NumRows);
            foreach (var wholeRow in data)
            {
                foreach (var colSign in wholeRow.Value)
                {
                    transpose.AddEntry(colSign.Key, wholeRow.Key, colSign.Value == 1);
                }
            }
            return transpose;
        }

        //TODO: do this in a dedicated Writer class
        public void WriteToConsole()
        {
            for (int i = 0; i < NumRows; ++i)
            {
                bool rowExists = data.TryGetValue(i, out Dictionary<int, int> colSigns);
                for (int j = 0; j < NumColumns; ++j)
                {
                    int val = 0;
                    if (rowExists) colSigns.TryGetValue(j, out val);
                    Console.Write($"{val,3}");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        //TODO: I think that it will pay off to transpose an all integer CSR matrix and store both. Especially in the case of 
        //      subdomain boolean matrices, that little extra memory should not be of concern.
        private Vector MultiplyTransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
            var result = new double[NumColumns];
            // Transpose it conceptually and multiply with the vector on the right. 
            foreach (var wholeRow in data)
            {
                foreach (var colSign in wholeRow.Value)
                {
                    result[colSign.Key] += colSign.Value * vector[wholeRow.Key];
                }
                
            }
            return Vector.CreateFromArray(result, false);
        }

        private Vector MultiplyUntransposed(Vector vector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
            var result = new double[NumRows];
            foreach (var wholeRow in data)
            {
                foreach (var colSign in wholeRow.Value)
                {
                    result[wholeRow.Key] += colSign.Value * vector[colSign.Key];
                }
            }
            return Vector.CreateFromArray(result, false);
        }

        
    }
}
