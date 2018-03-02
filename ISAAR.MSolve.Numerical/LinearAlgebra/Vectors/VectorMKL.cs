using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

//TODO: align data using mkl_malloc
//TODO: tensor product, vector2D, vector3D
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Vectors
{
    public class VectorMKL: IReducible
    {
        private readonly double[] data;

        private VectorMKL(double[] data)
        {
            this.data = data;
            this.Length = data.Length;
        }

        public int Length { get; }

        /// <summary>
        /// TODO: make this package-private. It should only be used for passing raw arrays to linear algebra libraries
        /// </summary>
        internal double[] InternalData { get { return data; } }

        public double this[int i]
        {
            get { return data[i]; }
            set { data[i] = value; }
        }

        public static VectorMKL CreateZero(int length)
        {
            return new VectorMKL(new double[length]);
        }

        public static VectorMKL CreateWithValue(int length, double value)
        {
            double[] data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = value;
            return new VectorMKL(data);
        }

        public static VectorMKL CreateFromArray(double[] data, bool copyArray = true)
        {
            if (copyArray)
            {
                double[] clone = new double[data.Length];
                Array.Copy(data, clone, data.Length);
                return new VectorMKL(clone);
            }
            else return new VectorMKL(data);
        }

        #region operators 
        public static VectorMKL operator +(VectorMKL v1, VectorMKL v2)
        {
            return v1.Axpy(1.0, v2);
        }

        public static VectorMKL operator -(VectorMKL v1, VectorMKL v2)
        {
            return v1.Axpy(-1.0, v2); //The order is important
        }

        public static VectorMKL operator *(double scalar, VectorMKL vector)
        {
            return vector.Scale(scalar);
        }

        public static double operator *(VectorMKL v1, VectorMKL v2)
        {
            return v1.DotProduct(v2); //Perhaps call BLAS directly
        }
        #endregion

        //Perhaps this should use BLAS
        public static VectorMKL CreateFromVector(VectorMKL vector)
        {
            double[] data = vector.data;
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new VectorMKL(clone);
        }

        /// <summary>
        /// result = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public VectorMKL Axpy(double scalar, VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        public void AxpyIntoThis(double scalar, VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            CBlas.Daxpy(Length, scalar, ref other.data[0], 1, ref this.data[0], 1);
        }

        //Perhaps this should use BLAS
        public double[] CopyToArray()
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return clone;
        }

        public double DotProduct(VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            return CBlas.Ddot(Length, ref this.data[0], 1, ref other.data[0], 1);
        }

        public bool Equals(VectorMKL other, ValueComparer comparer = null)
        {
            if (this.Length != other.Length) return false;
            if (comparer == null) comparer = new ValueComparer(1e-13);
            for (int i = 0; i < Length; ++i)
            {
                if (!comparer.AreEqual(this.data[i], other.data[i])) return false;
            }
            return true;
        }

        /// <summary>
        /// result = thisScalar * this + otherScalar * otherVector
        /// </summary>
        /// <returns></returns>
        public VectorMKL LinearCombination(double thisScalar, double otherScalar, VectorMKL otherVector)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Daxpby(Length, otherScalar, ref otherVector.data[0], 1, thisScalar, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationIntoThis(double thisScalar, double otherScalar, VectorMKL otherVector)
        {
            Preconditions.CheckVectorDimensions(this, otherVector);
            CBlas.Daxpby(Length, otherScalar, ref otherVector.data[0], 1, thisScalar, ref this.data[0], 1);
        }

        /// <summary>
        /// Computes the Hadamard product of two vectors: result[i] = this[i] * other[i]. This is not the dot product.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public VectorMKL MultiplyPointwise(VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = this.data[i] * other.data[i];
            return new VectorMKL(result);
        }

        /// <summary>
        /// Computes the Hadamard product of two vectors: this[i] = this[i] * other[i]. This is not the dot product.
        /// </summary>
        /// <param name="other"></param>
        public void MultiplyPointwiseIntoThis(VectorMKL other)
        {
            Preconditions.CheckVectorDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) this.data[i] *= other.data[i];
        }

        public double Norm2()
        {
            return CBlas.Dnrm2(Length, ref data[0], 1);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
        }

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public VectorMKL Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(Length, scalar, ref result[0], 1);
            return new VectorMKL(result);
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            CBlas.Dscal(Length, scalar, ref data[0], 1);
        }

        public void SetAll(double value)
        {
            for (int i = 0; i < Length; ++i) data[i] = value;
        }

        /// <summary>
        /// Returns a subvector containing only the entries at the provided indices
        /// </summary>
        /// <param name="indices">Indices of the entries to be returned. They must be 0 &lt; = i &lt; <see cref="Length"/>.</param>
        /// <returns></returns>
        public VectorMKL Slice(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new VectorMKL(subvector);
        }

        /// <summary>
        /// Returns a subvector containing the entries at the indices between the provided start (inclusive) and end (exclusive).
        /// </summary>
        /// <param name="startInclusive">The first index from which to copy entries.</param>
        /// <param name="endExclusive">The index after the last one until which to copy entries.</param>
        /// <returns></returns>
        public VectorMKL Slice(int startInclusive, int endExclusive)
        {
            int newLength = endExclusive - startInclusive;
            double[] subvector = new double[newLength];
            for (int i = 0; i < newLength; ++i) subvector[i] = data[startInclusive + i];
            return new VectorMKL(subvector);
        }

        public void WriteToConsole(Array1DFormatting format = null)
        {
            if (format == null) format = Array1DFormatting.Default;
            string separator = format.Separator;
            Console.Write(format.Start);
            for (int i = 0; i < Length; ++i)
            {
                Console.Write(separator + data[i]);
            }
            Console.WriteLine(separator + format.End);
        }

        /// <summary>
        /// Write the entries of the vector to a specified file. If the file doesn't exist a new one will be created.
        /// </summary>
        /// <param name="path">The path of the file and its extension.</param>
        /// <param name="append">If the file already exists: Pass <see cref="append"/> = true to write after the current end of 
        ///     the file. Pass<see cref="append"/> = false to overwrite the file.</param>
        /// <param name="format">Formatting options for how to print the vector entries.</param>
        public void WriteToFile(string path, bool append = false, Array1DFormatting format = null)
        {
            //TODO: incorporate this and WriteToConsole into a common function, where the user passes the stream and an object to 
            //deal with formating. Also add support for relative paths. Actually these methods belong in the "Logging" project, 
            // but since they are extremely useful they are implemented here for now.
            using (var writer = new StreamWriter(path, append))
            {
                string separator = format.Separator;
                writer.Write(format.Start);
                for (int i = 0; i < Length; ++i)
                {
                    writer.Write(separator + data[i]);
                }
                writer.WriteLine(separator + format.End);
            }
        }
    }
}
