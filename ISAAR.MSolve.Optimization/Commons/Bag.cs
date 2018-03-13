using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Commons
{
    /// <summary>
    /// Data structure with the following operations: 
    /// Constructor: O(n) running time, O(n) extra memory
    /// Get random entry: O(1) running time, no extra memory
    /// Remove specified entry: O(n) running time, O(n) extra memory
    /// </summary>
    /// <typeparam name="T">Can be anything</typeparam>
    class Bag<T>
    {
        private T[] entries;
        private readonly IGenerator rng;

        public Bag(IEnumerable<T> entries, IGenerator randomNumberGenerator = null)
        {
            this.entries = entries.ToArray();
            rng = (randomNumberGenerator != null) ? randomNumberGenerator : RandomNumberGenerationUtilities.troschuetzRandom;
        }

        public T GetRandom()
        {
            return entries[rng.Next(entries.Length)];
        }

        public bool Remove(T entry)
        {
            for (int i = 0; i < entries.Length; ++i)
            {
                if (entry.Equals(entries[i]))
                {
                    entries = CloneEntriesWithoutIndex(i);
                    return true;
                }
            }
            return false;
        }

        // Assumes the specified index is valid
        private T[] CloneEntriesWithoutIndex(int index)
        {
            T[] clone = new T[entries.Length - 1];
            Array.Copy(entries, clone, index);
            Array.Copy(entries, index + 1, clone, index, entries.Length - index - 1); // length = (N-1) - (index+1) + 1
            return clone;
        }
    }
}
