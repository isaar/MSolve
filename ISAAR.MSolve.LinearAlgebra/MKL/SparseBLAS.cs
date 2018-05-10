using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.MKL
{
    class SparseBLAS
    {
        //TODO: this has been deprecated. Does this hurt performance much? Should I use the inspector-executor functions instead?
        [DllImport("mkl_rt", CallingConvention = CallingConvention.Cdecl, EntryPoint = "mkl_dcsrcsc")]
        internal static extern int DCsrCsc(int[] job, int n, double[] acsr, int[] ja, int[] ia, 
            double[] acsc, int[] ja1, int[] ia1, ref int info);
    }
}
