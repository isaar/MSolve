using System;


namespace ISAAR.MSolve.Optimization.Commons
{
    public static class VectorOperations
    {
        public static double[] Add(double[] v1, double[] v2)
        {
            CheckDimensions(v1, v2);

            double[] result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i] + v2[i];
            }
            return result;
        }

        public static double[] Scale(double scalar, double[] vector)
        {
            double[] result = new double[vector.Length];

            for (int i = 0; i < vector.Length; i++)
            {
                result[i] = scalar * vector[i];
            }
            return result;
        }

        public static double[] Subtract(double[] v1, double[] v2)
        {
            CheckDimensions(v1, v2);

            double[] result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
            {
                result[i] = v1[i] - v2[i];
            }
            return result;
        }

        private static void CheckDimensions(double[] v1, double[] v2)
        {
            if (v1.Length != v2.Length)
            {
                throw new ArgumentException("Vectors must have equal lengths!");
            }
        }
    }
}
