using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Analyzers.Optimization.Commons
{
    /// <summary>
    /// A data structure that supports random removing of objects in O(1), but only after all items are inserted.
    /// Once the first remove is executed, insertion is locked and the entries are shuffled.
    /// It is not thread safe.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class ShuffledReadonlyBag<T>
    {
        private LinkedList<T> entries;
        private bool isShuffled;

        public ShuffledReadonlyBag()
        {
            this.entries = new LinkedList<T>();
            this.isShuffled = false;
        }

        public int Count { get { return entries.Count; } }

        public void Add(T item)
        {
            if (isShuffled) entries.AddLast(item);
            else throw new InvalidOperationException("Cannot add any more items after the first remove");
        }

        public T RemoveRandom()
        {
            if (!isShuffled)
            {
                Shuffle();
                isShuffled = true;
            }
            T item = entries.Last.Value;
            entries.RemoveLast();
            return item;
        }

        private void Shuffle()
        {
            //TODO: find a better way to implement the whole structure
            T[] asArray = new T[entries.Count];
            entries.CopyTo(asArray, 0);
            Random rng = new Random();
            int n = entries.Count;
            while (n > 1)
            {
                n--;
                int k = rng.Next(n + 1);
                T value = asArray[k];
                asArray[k] = asArray[n];
                asArray[n] = value;
            }
            entries = new LinkedList<T>(asArray);
        }
    }
}
