using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Utility methods to use when calling methods from <see cref="ILapackProvider"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class LapackUtilities
    {
        /// <summary>
        /// Use this value to initialize the output arguments named info. So far, negative values of the info are parameter 
        /// indices, thus the int.MinValue is safe to use as a default.  
        /// </summary>
        internal const int DefaultInfo = int.MinValue;

        /// <summary>
        /// Abstracts LAPACK subroutines where workspace array management is needed.
        /// </summary>
        internal delegate void RawLapackRoutine(double[] workArray, int offsetWork, int lengthWork);

        /// <summary>
        /// Most LAPACK functions return the same negative error codes for invalid parameters. Parameter <paramref name="info"/> 
        /// will not be checked to make sure it is negative.
        /// </summary>
        [Conditional("DEBUG")]
        internal static void ProcessNegativeInfo(int info)
        {
            if (info == DefaultInfo)
            {
                // first check the default info value, since it is negative.
                // info == default => the LAPACK call did not succeed. 
                throw new LapackException("Something went wrong with the LAPACK call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else
            {
                string suffix;
                if (info == -1) suffix = "st";
                else if (info == -2) suffix = "nd";
                else if (info == -3) suffix = "rd";
                else suffix = "th";
                throw new LapackException($"The {-info}{suffix} parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            // If positive values are possible they should be checked by the called
        }

        //TODO: The query info result is overwritten by the code from the second call and thus cannot be checked.
        /// <summary>
        /// For LAPACK subroutines where workspace array management is needed: First the subroutine is called to find a  
        /// reasonably good (for optimality ILAENV must be installed for each system) size for the workspace array. Then the 
        /// array is allocated and the subroutine is called again, but this time the actual operation is performed.
        /// </summary>
        internal static void QueryWorkspaceAndExecute(RawLapackRoutine routine)
        {
            // Tell LAPACK query to calculate the optimum workspace size, instead of executing the operation
            var workQuery = new double[1];
            int lengthWorkQuery = -1;
            routine(workQuery, 0, lengthWorkQuery);

            // Call LAPACK to actually perform the operation
            int lengthWorkOptim = (int)(workQuery[0]);
            if (lengthWorkOptim < 1) lengthWorkOptim = 1; //TODO: should I throw an exception instead?
            var workOptim = new double[lengthWorkOptim];
            routine(workOptim, 0, lengthWorkOptim);
        }
    }
}
