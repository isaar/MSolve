using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Length
{
    class ConstantIncrement2D: ICrackGrowthLengthLaw2D
    {
        private readonly double length;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="incrementLength">Usually 0.1 is used</param>
        public ConstantIncrement2D(double incrementLength)
        {
            if (incrementLength <= 0.0)
            {
                throw new ArgumentException("The increment must be positive, but was: " + incrementLength);
            }
            this.length = incrementLength;
        }

        public double ComputeGrowthLength(double sif1, double sif2)
        {
            return length;
        }
    }
}
