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
        /// result = vector1 + vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static VectorMKL Add(this VectorMKL vector1, VectorMKL vector2)
        {
            return vector1.Axpy(1.0, vector2);
        }

        /// <summary>
        /// vector1 = vector1 + vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        public static void AddIntoThis(this VectorMKL vector1, VectorMKL vector2)
        {
            vector1.AxpyIntoThis(1.0, vector2);
        }

        /// <summary>
        /// result = vector1 - vector1
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static VectorMKL Subtract(this VectorMKL vector1, VectorMKL vector2)
        {
            return vector1.Axpy(-1.0, vector2);
        }

        /// <summary>
        /// vector1 = vector1 - vector2
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        public static void SubtractIntoThis(this VectorMKL vector1, VectorMKL vector2)
        {
            vector1.AxpyIntoThis(-1.0, vector2);
        }

        public static VectorMKL[] RemoveDuplicatesFindMultiplicity()
        {
            throw new NotImplementedException("I kinda forgot why dimtsap wanted this...");
        }
    }
}
