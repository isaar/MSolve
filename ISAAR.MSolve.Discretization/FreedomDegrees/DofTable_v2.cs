
//TODO: This is for the case when we also number constrained dofs globally.

//using System.Collections.Generic;
//using ISAAR.MSolve.Numerical.Commons;

//namespace ISAAR.MSolve.Discretization.FreedomDegrees
//{
//    /// <summary>
//    /// A <see cref="ITable{TRow, TColumn, TValue}"/> that associates the freedom degrees of nodes with their ordinal number.
//    /// Authors: Serafeim Bakalakos
//    /// </summary>
//    /// <typeparam name="TDof">A freedom degree type.</typeparam>
//    public class DofTable_v2<TDof>: Table<IDiscretePoint, TDof, int> where TDof: IDof
//    {
//        public DofTable_v2(): base()
//        { }

//        private DofTable_v2(Dictionary<IDiscretePoint, Dictionary<TDof, int>> data): base(data)
//        { }

//        //TODO: this would be nice to have in Table too.
//        public DofTable_v2<TDof> DeepCopy()
//        {
//            var dataCopy = new Dictionary<IDiscretePoint, Dictionary<TDof, int>>();
//            foreach (var wholeRow in this.data)
//            {
//                // IDof and int are immutable, thus I can just copy the nested dictionary.
//                dataCopy.Add(wholeRow.Key, new Dictionary<TDof, int>(wholeRow.Value));
//            }
//            return new DofTable_v2<TDof>(dataCopy);
//        }

//        /// <summary>
//        /// Renumbers the dof indices according to the given permutation vector and direction. 
//        /// If (<paramref name="oldToNew"/> == true), then newIndex[dof] = <paramref name="permutation"/>[oldIndex[dof]].
//        /// Else oldIndex[dof] = <paramref name="permutation"/>[nwIndex[dof]]
//        /// </summary>
//        /// <param name="permutation">The permutation vector.</param>
//        /// <param name="oldToNew">The direction it should be applied to.</param>
//        public void Reorder(IReadOnlyList<int> permutation, bool oldToNew)
//        {
//            IReadOnlyList<int> permutationOldToNew;
//            if (oldToNew) permutationOldToNew = permutation;
//            else
//            {
//                var permutationArray = new int[permutation.Count];
//                for (int newIdx = 0; newIdx < permutation.Count; ++newIdx) permutationArray[permutation[newIdx]] = newIdx;
//                permutationOldToNew = permutationArray;
//            }

//            foreach (var nodeRow in data.Values)
//            {
//                var dofIDs = new List<KeyValuePair<TDof, int>>(nodeRow);
//                foreach (var dofIDPair in dofIDs)
//                {
//                    nodeRow[dofIDPair.Key] = permutationOldToNew[dofIDPair.Value];
//                }

//                // The following throws CollectionModified, although nothing dangerous happens 
//                //foreach (var dofID in nodeRow)
//                //{
//                //    nodeRow[dofID.Key] = permutationOldToNew[dofID.Value];
//                //}
//            }
//        }
//    }
//}
