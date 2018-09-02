using ISAAR.MSolve.LinearAlgebra.Exceptions;

// TODO: Perhaps I should have various info classes and wrap the MKL calls directly
namespace ISAAR.MSolve.LinearAlgebra.MKL
{
    /// <summary>
    /// Utility methods to use when calling Intel MKL functions.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class MklUtilities
    {
        /// <summary>
        /// Use this value to initialize the output arguments named info. So far, negative values of the info are parameter 
        /// indices, thus the int.MinValue is safe to use as a default.  
        /// </summary>
        public const int DefaultInfo = int.MinValue;

        /// <summary>
        /// Most MKL functions return the same negative error codes for invalid parameters. Parameter <paramref name="info"/> 
        /// will not be checked to make sure it is negative.
        /// </summary>
        /// <param name="info"></param>
        internal static MklException ProcessNegativeInfo(int info)
        {
            if (info == MklUtilities.DefaultInfo)
            {
                // first check the default info value, since it is negative.
                // info == default => the MKL call did not succeed. 
                return new MklException("Something went wrong with the MKL call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else
            {
                string suffix;
                if (info == -1) suffix = "st";
                else if (info == -2) suffix = "nd";
                else if (info == -3) suffix = "rd";
                else suffix = "th";
                return new MklException($"The {-info}{suffix} parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
        }
    }
}
