using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using System.Threading;
using System.IO;
using System.Globalization;
using System.Threading.Tasks;
using System.Collections.Specialized;
using System.Diagnostics;

namespace ISAAR.MSolve.Matrices
{
    public static class VectorExtensions
    {
        public static int AffinityCount = 0;

        public static void AssignTotalAffinityCount()
        {
            var b = new BitVector32((int)Process.GetCurrentProcess().ProcessorAffinity);
            int count = 0;
            for (int i = 0; i < 32; i++)
                if (b[1 << i]) count++;

            AffinityCount = count < 1 ? 1 : count;
        }

        public static void AssignTotalAffinityCount(int count)
        {
            AffinityCount = count < 1 ? 1 : count;
        }

        private static IEnumerable<Tuple<int, int, int>> GetVectorLimits(int size, int chunks)
        {
            int chunkSize = ((size % chunks) == 0) ? size / chunks : ((int)((size) / ((float)chunks)) + 1);
            int currentChunk = 0;
            int endPos = 0;
            while (currentChunk < chunks)
            {
                int chunk = Math.Min(chunkSize, size - currentChunk * chunkSize);
                endPos += chunk;
                currentChunk++;
                yield return new Tuple<int, int, int>(currentChunk - 1, endPos - chunk, endPos);
            }
        }

        public static IEnumerable<Tuple<int, int, int>> PartitionLimits<T>(this T[] vector, int chunks)
        {
            int size = vector.Length;
            return GetVectorLimits(size, chunks);
        }

        public static IEnumerable<Tuple<int, int, int>> PartitionLimits(this double[] vector, int chunks)
        {
            int size = vector.Length;
            return GetVectorLimits(size, chunks);
        }

        public static IEnumerable<Tuple<int, int, int>> PartitionLimits(this Vector<double> vector, int chunks)
        {
            int size = vector.Length;
            return GetVectorLimits(size, chunks);
        }
    }

    public class Vector<T> : IVector<T>
    {
        private readonly T[] data;

        public Vector(int length)
        {
            data = new T[length];
        }

        public Vector(T[] data)
        {
            this.data = data;
        }

        public T[] Data
        {
            get { return data; }
        }

        public double Norm
        {
            get
            {
                if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
                double norm = 0;
                double[] v1Data = data as double[];
                for (int i = 0; i < this.Length; i++) norm += v1Data[i] * v1Data[i];
                return Math.Sqrt(norm);
            }
        }

        public static Vector<T> operator *(double s1, Vector<T> v2)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] v2Data = v2.data as double[];
            double[] result = new double[v2.Length];
            for (int i = 0; i < v2.Length; i++) result[i] = s1 * v2Data[i];
            return (new Vector<double>(result)) as Vector<T>;
        }

        private static double DoDot(int segment, double[] d1, double[] d2)
        {
            int l = (int)(d1.Length / VectorExtensions.AffinityCount);
            int start = segment * l;
            int finish = (segment + 1) * l;
            if (segment == VectorExtensions.AffinityCount - 1) finish = d1.Length;
            double result = 0;
            for (int i = start; i < finish; i++) result += d1[i] * d2[i];
            return result;
        }

        public static double operator *(Vector<T> v1, Vector<T> v2)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double result = 0;
            double[] v1Data = v1.data as double[];
            double[] v2Data = v2.data as double[];
            //for (int i = 0; i < v1.Length; i++) result += v1Data[i] * v2Data[i];

            int iProcs = VectorExtensions.AffinityCount;
            double[] results = new double[iProcs];
            Parallel.For(0, iProcs, i => { results[i] = DoDot(i, v1Data, v2Data); });
            for (int i = 0; i < iProcs; i++) result += results[i];
            return result;
        }

        public static T[] operator +(Vector<T> v1, Vector<T> v2)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] result = new double[v1.Length];
            double[] v1Data = v1.data as double[];
            double[] v2Data = v2.data as double[];
            for (int i = 0; i < v1.Length; i++) result[i] = v1Data[i] + v2Data[i];
            return result as T[];
        }

        public static T[] operator -(Vector<T> v1, Vector<T> v2)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] result = new double[v1.Length];
            double[] v1Data = v1.data as double[];
            double[] v2Data = v2.data as double[];
            for (int i = 0; i < v1.Length; i++) result[i] = v1Data[i] - v2Data[i];
            return result as T[];
        }

        public void Scale(double scale)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] vData = data as double[];
            for (int i = 0; i < data.Length; i++) vData[i] *= scale;
        }

        public void Add(Vector<T> v)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] v1Data = data as double[];
            double[] v2Data = v.Data as double[];
            for (int i = 0; i < data.Length; i++) v1Data[i] += v2Data[i];
        }

        public void Subtract(Vector<T> v)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] v1Data = data as double[];
            double[] v2Data = v.Data as double[];
            for (int i = 0; i < data.Length; i++) v1Data[i] -= v2Data[i];
        }
        #region IVector<T> Members

        public int Length
        {
            get { return data.Length; }
        }

        public T this[int x]
        {
            get { return data[x]; }
            set { data[x] = value; }
        }

        public void Multiply(double coefficient)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] vData = data as double[];
            for (int i = 0; i < this.Length; i++) vData[i] *= (double)coefficient;
        }

        public void CopyTo(Array array, int index)
        {
            data.CopyTo(array, index);
        }

        public void Clear()
        {
            Array.Clear(data, 0, data.Length);
        }

        public void WriteToFile(string name)
        {
            if (typeof(T) != typeof(double)) throw new InvalidOperationException("Only double type is supported.");
            double[] mData = data as double[];

            StreamWriter sw = new StreamWriter(name);
            foreach (double d in mData)
                sw.WriteLine(d.ToString("g17", new CultureInfo("en-US", false).NumberFormat));
            sw.Close();
        }

        public Vector<double>[] RemoveDuplicatesFindMultiplicity()
        {
            throw new NotImplementedException();
        }

        public Vector<double> FindUnionWithVector(Vector<double> vector)
        {
            throw new NotImplementedException();
        }

        public Vector<double> FindIntersectionWithVector(Vector<double> vector)
        {
            throw new NotImplementedException();
        }

        public void SortAscending()
        {
            Array.Sort(data);
        }

        public void SortDescending()
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
