using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    class DOFTable<TDOF>: Table<XNode2D, TDOF, int> where TDOF: IDOF
    {
        public void Reorder(IReadOnlyList<int> permutationOldToNew)
        {
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
