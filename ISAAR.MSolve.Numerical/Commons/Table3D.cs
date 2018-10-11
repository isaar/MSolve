using System.Collections;
using System.Collections.Generic;

namespace ISAAR.MSolve.Numerical.Commons
{
    /// <summary>
    /// Interface for table data structures, which associate ordered pairs (row, column) with values. Contrary to matrices, not 
    /// all (row, column) pairs of a table need to be associated with values.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TRow"></typeparam>
    /// <typeparam name="TColumn"></typeparam>
    /// <typeparam name="TValue"></typeparam>
    public class Table3D<TKey1, TKey2, TKey3, TValue> : ITable3D<TKey1, TKey2, TKey3, TValue>
    {
        private const int defaultInitialCapacity = 1; //There will be at least 1. TODO: perhaps get this from Dictionary class.
        private readonly int initialCapacityForEachDim;
        private readonly Dictionary<TKey1, Dictionary<TKey2, Dictionary<TKey3, TValue>>> data;

        public Table3D(int initialCapacityForEachDim = defaultInitialCapacity)
        {
            this.initialCapacityForEachDim = initialCapacityForEachDim;
            this.data = new Dictionary<TKey1, Dictionary<TKey2, Dictionary<TKey3, TValue>>>(initialCapacityForEachDim);
        }

        public int EntryCount//TODO: perhaps this should be cached somehow
        {
            get
            {
                int count = 0;
                foreach (var key1Entries in data)
                {
                    foreach (var key2Entries in key1Entries.Value) count += key2Entries.Value.Count;
                }
                return count;
            }
        }

        public TValue this[TKey1 key1, TKey2 key2, TKey3 key3]
        {
            get => data[key1][key2][key3];

            set
            {
                bool containsKey1 = data.TryGetValue(key1, out Dictionary<TKey2, Dictionary<TKey3, TValue>> key1Entries);
                if (!containsKey1)
                {
                    key1Entries = new Dictionary<TKey2, Dictionary<TKey3, TValue>>(initialCapacityForEachDim);
                    data.Add(key1, key1Entries);
                }

                bool containsKey2 = key1Entries.TryGetValue(key2, out Dictionary<TKey3, TValue> key2Entries);
                if (!containsKey2)
                {
                    key2Entries = new Dictionary<TKey3, TValue>(initialCapacityForEachDim);
                    key1Entries.Add(key2, key2Entries);
                }

                key2Entries[key3] = value; // This allows changing the value after an entry has been added.
            }
        }

        public bool Contains(TKey1 key1, TKey2 key2, TKey3 key3)
        {
            bool containsKey1 = data.TryGetValue(key1, out Dictionary<TKey2, Dictionary<TKey3, TValue>> key1Entries);
            if (!containsKey1) return false;

            bool containsKey2 = key1Entries.TryGetValue(key2, out Dictionary<TKey3, TValue> key2Entries);
            if (!containsKey2) return false;

            return key2Entries.ContainsKey(key3);
        }

        public IEnumerator<(TKey1 key1, TKey2 key2, TKey3 key3, TValue val)> GetEnumerator()
        {
            foreach (var key1Entries in data)
            {
                foreach (var key2Entries in key1Entries.Value)
                {
                    foreach (var key3Entries in key2Entries.Value)
                    {
                        yield return (key1Entries.Key, key2Entries.Key, key3Entries.Key, key3Entries.Value);
                    }
                }
            }
        }

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        public bool TryGetValue(TKey1 key1, TKey2 key2, TKey3 key3, out TValue value)
        {
            bool containsKey1 = data.TryGetValue(key1, out Dictionary<TKey2, Dictionary<TKey3, TValue>> key1Entries);
            if (!containsKey1)
            {
                value = default(TValue);
                return false;
            }

            bool containsKey2 = key1Entries.TryGetValue(key2, out Dictionary<TKey3, TValue> key2Entries);
            if (!containsKey2)
            {
                value = default(TValue);
                return false;
            }

            return key2Entries.TryGetValue(key3, out value);
        }
    }
}
