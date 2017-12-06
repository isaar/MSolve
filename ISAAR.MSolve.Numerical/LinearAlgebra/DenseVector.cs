using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    class DenseVector: IVector
    {
        private readonly double[] data;

        public double this[int index]
        {
            get { return data[index]; }
            set { data[index] = value; }
        }

        public int Length { get { return data.Length; } }

        #region construction
        public static DenseVector CreateZero(int length)
        {
            return new DenseVector(new double[length]);
        }

        public static DenseVector CreateWithValue(int length, double value)
        {
            double[] data = new double[length];
            for (int i = 0; i < length; ++i) data[i] = value;
            return new DenseVector(data);
        }

        public static DenseVector CreateFromArray(double[] data, bool copyArray = true)
        {
            if (copyArray)
            {
                double[] clone = new double[data.Length];
                for (int i = 0; i < data.Length; ++i) clone[i] = data[i];
                return new DenseVector(clone);
            }
            else return new DenseVector(data);
        }

        public static DenseVector CreateFromVector(IVectorView vector)
        {
            return new DenseVector(vector.CopyToArray());
        }

        private DenseVector(double[] data)
        {
            this.data = data;
        }
        #endregion

        public IVector Add(IVectorView vector)
        {
            return DoPointwise(vector, (x1, x2) => x1 + x2);
        }

        public void AddIntoThis(IVectorView vector)
        {
            DoPointwiseIntoThis(vector, (x1, x2) => x1 + x2);
        }

        public double[] CopyToArray()
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return clone;
        }

        public IVector DoPointwise(IVectorView vector, Func<double, double, double> operation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = operation(data[i], vector[i]);
            return new DenseVector(result);
        }

        public void DoPointwiseIntoThis(IVectorView vector, Func<double, double, double> operation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            for (int i = 0; i < data.Length; ++i) data[i] = operation(data[i], vector[i]);
        }

        public IVector DoToAllEntries(Func<double, double> operation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = operation(data[i]);
            return new DenseVector(result);
        }

        public IVector ExtractSubvector(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new DenseVector(subvector);
        }

        public void DoToAllEntriesIntoThis(Func<double, double> operation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = operation(data[i]);
        }

        public double MultiplyDot(IVectorView vector)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double result = 0.0;
            for (int i = 0; i < data.Length; ++i) result += data[i] * vector[i];
            return result;
        }

        public IVector MultiplyPointwise(IVectorView vector)
        {
            return DoPointwise(vector, (x1, x2) => x1 * x2);
        }

        public void MultiplyPointwiseIntoThis(IVectorView vector)
        {
            DoPointwiseIntoThis(vector, (x1, x2) => x1 * x2);
        }

        public IVector MultiplyScalar(double scalar)
        {
            return DoToAllEntries(x => x * scalar);
        }

        public void MultiplyScalarIntoThis(double scalar)
        {
            DoToAllEntriesIntoThis(x => x * scalar);
        }

        public void Print()
        {
            for (int i = 0; i < data.Length-2; ++i)
            {
                Console.Write(data[i]);
                Console.Write(' ');
            }
            Console.Write(data[data.Length-1]);
        }

        DenseVector[] RemoveDuplicatesFindMultiplicity()
        {
            throw new NotImplementedException();
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, 
            ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
        }

        public void SetAll(double value)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = value; // Awkward with DoToAllEntriesIntoThis
        }

        public IVector Subtract(IVectorView vector)
        {
            return DoPointwise(vector, (x1, x2) => x1 - x2);
        }

        public void SubtractIntoThis(IVectorView vector)
        {
            DoPointwiseIntoThis(vector, (x1, x2) => x1 - x2);
        }

        public void Write(string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < data.Length - 2; ++i)
            {
                writer.Write(data[i]);
                writer.Write(' ');
            }
            writer.Write(data[data.Length - 1]);
        }
    }
}
