using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Utility methods for parallel vector operations.
    /// Authors: George Stavroulakis
    /// </summary>
    public static class ParallelExtensions
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

        public static IEnumerable<Tuple<int, int, int>> PartitionLimits(this IVectorView vector, int chunks)
        {
            int size = vector.Length;
            return GetVectorLimits(size, chunks);
        }
    }
}
