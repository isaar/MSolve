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

        public double[] CopyToArray()
        {
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return clone;
        }

        public IVector DoPointwise(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(data[i], vector[i]);
            return new DenseVector(result);
        }

        public void DoPointwiseIntoThis(IVectorView vector, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            for (int i = 0; i < data.Length; ++i) data[i] = binaryOperation(data[i], vector[i]);
        }

        public IVector DoToAllEntries(Func<double, double> unaryOperation)
        {
            double[] result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = unaryOperation(data[i]);
            return new DenseVector(result);
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
        }

        public IVector ExtractSubvector(int[] indices)
        {
            double[] subvector = new double[indices.Length];
            for (int i = 0; i < indices.Length; ++i) subvector[i] = data[indices[i]];
            return new DenseVector(subvector);
        }

        public double MultiplyDot(IVectorView vector)
        {
            Preconditions.CheckVectorDimensions(this, vector);
            double result = 0.0;
            for (int i = 0; i < data.Length; ++i) result += data[i] * vector[i];
            return result;
        }

        public void Print()
        {
            for (int i = 0; i < data.Length; ++i) 
            {
                Console.WriteLine(data[i]);
            }
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
            for (int i = 0; i < data.Length; ++i) data[i] = value; 
        }

        public void Write(string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < data.Length; ++i)
            {
                writer.WriteLine(data[i]);
            }
        }
    }
}
