using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//TODO: Move the operators here when C# supports extension operators
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
{
    public static class VectorExtensionsMKL
    {
        /// <summary>
        /// result = v1 + v2
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static VectorMKL Add(this VectorMKL v1, VectorMKL v2)
        {
            return v1.Axpy(1.0, v2);
        }

        /// <summary>
        /// v1 = v1 + v2
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        public static void AddIntoThis(this VectorMKL v1, VectorMKL v2)
        {
            v1.AxpyIntoThis(1.0, v2);
        }

        /// <summary>
        /// result = v1 - v2
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static VectorMKL Subtract(this VectorMKL v1, VectorMKL v2)
        {
            return v1.Axpy(-1.0, v2);
        }

        /// <summary>
        /// v1 = v1 - v2
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        public static void SubtractIntoThis(this VectorMKL v1, VectorMKL v2)
        {
            v1.AxpyIntoThis(-1.0, v2);
        }

        public static VectorMKL[] RemoveDuplicatesFindMultiplicity()
        {
            throw new NotImplementedException("I kinda forgot why dimtsap wanted this...");
        }
    }
}
