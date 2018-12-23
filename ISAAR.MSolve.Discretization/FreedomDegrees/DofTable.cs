using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// A <see cref="ITable{TRow, TColumn, TValue}"/> that associates the freedom degrees of nodes with their ordinal number.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TDof">A freedom degree type.</typeparam>
    public class DofTable: Table<INode, DOFType, int>
    {
        public DofTable(): base()
        { }

        private DofTable(Dictionary<INode, Dictionary<DOFType, int>> data): base(data)
        { }

        //TODO: this would be nice to have in Table too.
        public DofTable DeepCopy()
        {
            var dataCopy = new Dictionary<INode, Dictionary<DOFType, int>>();
            foreach (var wholeRow in this.data)
            {
                // IDof and int are immutable, thus I can just copy the nested dictionary.
                dataCopy.Add(wholeRow.Key, new Dictionary<DOFType, int>(wholeRow.Value));
            }
            return new DofTable(dataCopy);
        }

        /// <summary>
        /// Renumbers the dof indices according to the given permutation vector and direction. 
        /// If (<paramref name="oldToNew"/> == true), then newIndex[dof] = <paramref name="permutation"/>[oldIndex[dof]].
        /// Else oldIndex[dof] = <paramref name="permutation"/>[nwIndex[dof]]
        /// </summary>
        /// <param name="permutation">The permutation vector.</param>
        /// <param name="oldToNew">The direction it should be applied to.</param>
        internal void Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            IReadOnlyList<int> permutationOldToNew;
            if (oldToNew) permutationOldToNew = permutation;
            else
            {
                //TODO: is it necessary to create a temp oldToNew array?
                var permutationArray = new int[permutation.Count];
                for (int newIdx = 0; newIdx < permutation.Count; ++newIdx) permutationArray[permutation[newIdx]] = newIdx;
                permutationOldToNew = permutationArray;
            }

            foreach (var nodeRow in data.Values)
            {
                var dofIDs = new List<KeyValuePair<DOFType, int>>(nodeRow);
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

        //TODO: use Table.ModifyValues() for this.
        internal void ReorderNodeMajor(IReadOnlyList<INode> sortedNodes)
        {
            int dofIdx = -1;
            foreach (INode node in sortedNodes)
            {
                bool isNodeContained = data.TryGetValue(node, out Dictionary<DOFType, int> dofsOfNode);
                if (isNodeContained)
                {
                    // We cannot not update the dictionary during iterating it, thus we copy its Key collection to a list first.
                    foreach (DOFType dofType in dofsOfNode.Keys.ToList()) dofsOfNode[dofType] = ++dofIdx;
                }
            }
        }

        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            foreach (var nodeData in data.OrderBy(entry => entry.Key.ID))
            {
                builder.AppendLine($"Node {nodeData.Key.ID}:");
                foreach (var dofPair in nodeData.Value)
                {
                    builder.Append("\t");
                    builder.AppendLine($"Dof type = {dofPair.Key} - Global index = {dofPair.Value}");
                }
            }
            return builder.ToString();
        }
    }
}
