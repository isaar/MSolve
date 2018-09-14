using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Tests
{
    internal static class Utilities
    {
        internal static bool AreMatricesEqual(IMatrix2D matrix1, IMatrix2D matrix2, double tolerance)
        {
            if ((matrix1.Rows != matrix2.Rows) || (matrix1.Columns != matrix2.Columns)) return false;
            for (int i = 0; i < matrix1.Rows; ++i)
            {
                for (int j = 0; j < matrix1.Columns; ++j)
                {
                    if (!AreValuesEqual(matrix1[i, j], matrix2[i, j], tolerance)) return false;
                }
            }
            return true;
        }

        internal static bool AreTensorsEqual(IReadOnlyList<double[]> tensors1, IReadOnlyList<double[]> tensors2, double tolerance)
        {
            if (tensors1.Count != tensors2.Count) return false;
            for (int i = 0; i < tensors1.Count; ++i)
            {
                if (tensors1[i].Length != tensors2[i].Length) return false;
                for (int j = 0; j < tensors1[i].Length; ++j)
                {
                    if (!AreValuesEqual(tensors1[i][j], tensors2[i][j], tolerance)) return false;
                }
            }
            return true;
        }

        internal static bool AreValuesEqual(double value1, double value2, double tolerance)
        {
            if (Math.Abs(value2) <= tolerance) // Can't divide with expected ~= 0. 
            {
                if (Math.Abs(value1) <= tolerance) return true;
                else return false;
            }
            else return (Math.Abs(1.0 - value1 / value2) < tolerance) ? true : false;
        }
    }
}
