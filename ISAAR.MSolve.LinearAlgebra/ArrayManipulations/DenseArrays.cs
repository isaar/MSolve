using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.ArrayManipulations
{
    static class DenseArrays
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="tolerance">Can be zero</param>
        /// <returns></returns>
        public static bool IsZero(double[] array, double tolerance)
        {
            if (tolerance == 0)
            {
                for (int i = 0; i < array.Length; ++i)
                {
                    if (array[i] != 0.0) return false;
                }
                return true;
            }
            else
            {
                for (int i = 0; i < array.Length; ++i)
                {
                    if (Math.Abs(array[i]) > tolerance) return false;
                }
                return true;
            }
        }
    }
}
