using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Commons
{
    public static class CollectionExtensions
    {
        /// <summary>
        /// Searches for the specified object and returns the zero-based index of the first occurrence within the entire list.
        /// If none is found, -1 is returned.
        /// </summary>
        /// <typeparam name="TValue"></typeparam>
        /// <param name="list"></param>
        /// <param name="value">The object to search for.</param>
        /// <remarks>
        /// The list is searched forward starting at the first element and ending at the last element. This method determines 
        /// equality using <see cref="object.Equals(object)"/>. This method performs a linear search; therefore, this 
        /// method is an O(n) operation, where n is Count.
        /// </remarks>
        public static int FindFirstIndex<TValue>(this IReadOnlyList<TValue> list, TValue value)
        {
            for (int i = 0; i < list.Count; ++i)
            {
                if (list[i].Equals(value)) return i;
            }
            return -1;
        }

        /// <summary>
        /// Searches for an element of the list that satisfies a given condition and returns the zero-based index of the first 
        /// occurrence within the entire list. If none is found, -1 is returned.
        /// </summary>
        /// <typeparam name="TValue"></typeparam>
        /// <param name="list"></param>
        /// <param name="value">The object to search for.</param>
        /// <param name="predicate">A function to test each element for a condition.</param>
        /// <remarks>
        /// The list is searched forward starting at the first element and ending at the last element. This method performs a 
        /// linear search; therefore, this method is an O(n) operation, where n is Count.
        /// </remarks>
        public static int FindFirstIndex<TValue>(this IReadOnlyList<TValue> list, TValue value, Func<TValue, bool> predicate)
        {
            for (int i = 0; i < list.Count; ++i)
            {
                if (predicate(list[i])) return i;
            }
            return -1;
        }

        /// <summary>
        /// Searches for the specified object and returns the zero-based index of the last occurrence within the entire list.
        /// If none is found, -1 is returned.
        /// </summary>
        /// <typeparam name="TValue"></typeparam>
        /// <param name="list"></param>
        /// <param name="value">The object to search for.</param>
        /// <remarks>
        /// The list is searched backward starting at the last element and ending at the first element. This method determines 
        /// equality using <see cref="object.Equals(object)"/>. This method performs a linear search; therefore, this 
        /// method is an O(n) operation, where n is Count.
        /// </remarks>
        public static int FindLastIndex<TValue>(this IReadOnlyList<TValue> list, TValue value)
        {
            for (int i = list.Count - 1; i >= 0; --i)
            {
                if (list[i].Equals(value)) return i;
            }
            return -1;
        }

        /// <summary>
        /// Searches for an element of the list that satisfies a given condition and returns the zero-based index of the last 
        /// occurrence within the entire list. If none is found, -1 is returned.
        /// </summary>
        /// <typeparam name="TValue"></typeparam>
        /// <param name="list"></param>
        /// <param name="value">The object to search for.</param>
        /// <param name="predicate">A function to test each element for a condition.</param>
        /// <remarks>
        /// The list is searched backward starting at the last element and ending at the first element. This method performs a 
        /// linear search; therefore, this method is an O(n) operation, where n is Count.
        /// </remarks>
        public static int FindLastIndex<TValue>(this IReadOnlyList<TValue> list, TValue value, Func<TValue, bool> predicate)
        {
            for (int i = list.Count - 1; i >= 0; --i)
            {
                if (predicate(list[i])) return i;
            }
            return -1;
        }
    }
}
