using ISAAR.MSolve.LinearAlgebra.Exceptions;

// TODO: Perhaps I should have various info classes and wrap the LAPACK calls directly
namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    /// <summary>
    /// Utility methods to use when calling LAPACK functions.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class LapackUtilities
    {
        /// <summary>
        /// Use this value to initialize the output arguments named info. So far, negative values of the info are parameter 
        /// indices, thus the int.MinValue is safe to use as a default.  
        /// </summary>
        public const int DefaultInfo = int.MinValue;

        /// <summary>
        /// Most LAPACK functions return the same negative error codes for invalid parameters. Parameter <paramref name="info"/> 
        /// will not be checked to make sure it is negative.
        /// </summary>
        /// <param name="info"></param>
        internal static LapackException ProcessNegativeInfo(int info)
        {
            if (info == LapackUtilities.DefaultInfo)
            {
                // first check the default info value, since it is negative.
                // info == default => the LAPACK call did not succeed. 
                return new LapackException("Something went wrong with the LAPACK call."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
            else
            {
                string suffix;
                if (info == -1) suffix = "st";
                else if (info == -2) suffix = "nd";
                else if (info == -3) suffix = "rd";
                else suffix = "th";
                return new LapackException($"The {-info}{suffix} parameter has an illegal value."
                    + " Please contact the developer responsible for the linear algebra project.");
            }
        }
    }
}
