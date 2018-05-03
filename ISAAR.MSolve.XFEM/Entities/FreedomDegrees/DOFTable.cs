using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    class DOFTable<TDOF>: Table<XNode2D, TDOF, int> where TDOF: IDOF
    {
        public DOFTable(): base()
        { }

        private DOFTable(Dictionary<XNode2D, Dictionary<TDOF, int>> data): base(data)
        { }

        //TODO: this would be nice to have in Table too.
        public DOFTable<TDOF> DeepCopy()
        {
            var dataCopy = new Dictionary<XNode2D, Dictionary<TDOF, int>>();
            foreach (var wholeRow in this.data)
            {
                // IDOF and int are immutable, thus I can just copy the nested dictionary.
                dataCopy.Add(wholeRow.Key, new Dictionary<TDOF, int>(wholeRow.Value));
            }
            return new DOFTable<TDOF>(dataCopy);
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
                var dofIDs = new List<KeyValuePair<TDOF, int>>(nodeRow);
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
