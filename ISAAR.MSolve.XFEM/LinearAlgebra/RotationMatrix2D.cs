using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.LinearAlgebra
{
    // TODO: The matrix-matrix and matrix-vector operations should be extracted into a generic rectangular dense 
    // matrix. Rotation matrix classes should either be removed or just handle their construction (this could be the 
    // responsibility of a static class). Alternatively, I could have them as explicit IMatrix implementations and have 
    // optimized operations (e.g unrolled loops, cache friendly tweaks)
    class RotationMatrix2D
    {
        private readonly double[,] data = new double[2,2];

        public double this[int row, int col]
        {
            get { return data[row, col]; }
        }

        public RotationMatrix2D(double rotationAngle)
        {
            data[0, 0] = Math.Cos(rotationAngle);
            data[0, 1] = Math.Sin(rotationAngle);
            data[1, 0] = -data[0, 1];
            data[1, 1] = data[0, 0];
        }

        private RotationMatrix2D(double[,] data)
        {
            this.data = data;
        }

        public RotationMatrix2D Transpose()
        {
            double[,] transpose = new double[2,2];
            // Unrolled for more efficiency
            transpose[0, 0] = data[0, 0];
            transpose[0, 1] = data[1, 0];
            transpose[1, 0] = data[0, 1];
            transpose[1, 1] = data[1, 1];
            return new RotationMatrix2D(transpose);
        }

        public RotationMatrix2D Inverse()
        {
            return Transpose();
        }

        /// <summary>
        /// Performs the following matrix-vector multiplication: this * columnVector. Doesn't check dimensions.
        /// </summary>
        /// <param name="columnVector">A 2x1 column vector</param>
        /// <returns></returns>
        public double[] MultiplyRight(double[] columnVector)
        {
            double[] result = new double[2];
            // Unrolled for more efficiency
            result[0] = data[0, 0] * columnVector[0] + data[0, 1] * columnVector[1];
            result[1] = data[1, 0] * columnVector[0] + data[1, 1] * columnVector[1];
            return result;
        }

        /// <summary>
        /// Performs the following matrix multiplication: this * matrix. Doesn't check dimensions.
        /// </summary>
        /// <param name="matrix">A 2x2 matrix</param>
        /// <returns></returns>
        public double[,] MultiplyRight(double[,] matrix)
        {
            double[,] result = new double[2,2];
            // Unrolled for more efficiency
            result[0, 0] = data[0, 0] * matrix[0, 0] + data[0, 1] * matrix[1, 0];
            result[0, 1] = data[0, 0] * matrix[0, 1] + data[0, 1] * matrix[1, 1];
            result[1, 0] = data[1, 0] * matrix[0, 0] + data[1, 1] * matrix[1, 0];
            result[1, 1] = data[1, 0] * matrix[0, 1] + data[1, 1] * matrix[1, 1];
            return result;
        }

        /// <summary>
        /// Performs the following matrix-vector multiplication: rowVector * this. Doesn't check dimensions.
        /// </summary>
        /// <param name="matrix">A 1x2 row vector</param>
        /// <returns></returns>
        public double[] MultiplyLeft(double[] rowVector)
        {
            double[] result = new double[2];
            // Unrolled for more efficiency
            result[0] = rowVector[0] * data[0, 0] + rowVector[1] * data[1, 0];
            result[1] = rowVector[0] * data[0, 1] + rowVector[1] * data[1, 1];
            return result;
        }

        /// <summary>
        /// Performs the following matrix multiplication: matrix * this. Doesn't check dimensions.
        /// </summary>
        /// <param name="matrix">A 2x2 matrix</param>
        /// <returns></returns>
        public double[,] MultiplyLeft(double[,] matrix)
        {
            double[,] result = new double[2, 2];
            // Unrolled for more efficiency
            result[0, 0] = matrix[0, 0] * data[0, 0] + matrix[0, 1] * data[1, 0];
            result[0, 1] = matrix[0, 0] * data[0, 1] + matrix[0, 1] * data[1, 1];
            result[1, 0] = matrix[1, 0] * data[0, 0] + matrix[1, 1] * data[1, 0];
            result[1, 1] = matrix[1, 0] * data[0, 1] + matrix[1, 1] * data[1, 1];
            return result;
        }
    }
}
