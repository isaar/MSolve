using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;

namespace ISAAR.MSolve.LinearAlgebra.Reordering
{
    public class OrderingAMD
    {
        public OrderingAMD()
        {
            //TODO: add options for dense rows and aggressive absorption.
        }

        public int[] FindPermutation(int order, int nonZerosUpper, int[] cscRowIndices, int[] cscColOffsets)
        {
            var permutation = new int[order];
            IntPtr common = SuiteSparseUtilities.CreateCommon(0, 0);
            if (common == IntPtr.Zero) throw new SuiteSparseException("Failed to initialize SuiteSparse.");
            int status = SuiteSparseUtilities.ReorderAMDUpper(order, nonZerosUpper, cscRowIndices, cscColOffsets, permutation, 
                common);
            if (status == 0) throw new SuiteSparseException("AMD failed. This could be caused by the matrix being so large it"
                + " cannot be processed with the available memory.");
            SuiteSparseUtilities.DestroyCommon(ref common);
            return permutation;
        }
    }
}
