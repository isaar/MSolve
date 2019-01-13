using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.FreedomDegrees.Ordering
{
    class DofTable<TDof>: Table<XNode2D, TDof, int> where TDof: IDof
    {
        public DofTable(): base()
        { }

        private DofTable(Dictionary<XNode2D, Dictionary<TDof, int>> data): base(data)
        { }

        //TODO: this would be nice to have in Table too.
        public DofTable<TDof> DeepCopy()
        {
            var dataCopy = new Dictionary<XNode2D, Dictionary<TDof, int>>();
            foreach (var wholeRow in this.data)
            {
                // IDof and int are immutable, thus I can just copy the nested dictionary.
                dataCopy.Add(wholeRow.Key, new Dictionary<TDof, int>(wholeRow.Value));
            }
            return new DofTable<TDof>(dataCopy);
        }

        /// <summary>
        /// Renumbers the dof indices according th the given permutation vector and direction. 
        /// If (<paramref name="oldToNew"/> == true), then newIndex[dof] = <paramref name="permutation"/>[oldIndex[dof]].
        /// Else oldIndex[dof] = <paramref name="permutation"/>[nwIndex[dof]]
        /// </summary>
        /// <param name="permutation">The permutation vector.</param>
        /// <param name="oldToNew">The direction it should be applied to.</param>
        public void Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            IReadOnlyList<int> permutationOldToNew;
            if (oldToNew) permutationOldToNew = permutation;
            else
            {
                var permutationArray = new int[permutation.Count];
                for (int newIdx = 0; newIdx < permutation.Count; ++newIdx) permutationArray[permutation[newIdx]] = newIdx;
                permutationOldToNew = permutationArray;
            }

            foreach (var nodeRow in data.Values)
            {
                var dofIDs = new List<KeyValuePair<TDof, int>>(nodeRow);
                foreach (var dofIDPair in dofIDs)
                {
                    nodeRow[dofIDPair.Key] = permutationOldToNew[dofIDPair.Value];
                }

                // The following throws CollectionModified, although nothing dangerous happens 
                //foreach (var dofID in nodeRow)
                //{
                //    nodeRow[dofID.Key] = permutationOldToNew[dofID.Value];
                //}
            }
        }
    }
}
